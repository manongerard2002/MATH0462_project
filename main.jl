using Pkg
Pkg.activate(".")
Pkg.instantiate() #only required the first time you execute the code
using MATH0462_project, MATH0462_project.Reversibles
using JuMP, Gurobi, Printf, Random, JSON3

Random.seed!(42)

MINUTE = 10
SOFT_TIMEOUT = (MINUTE - 1) * 60
HARD_TIMEOUT = (MINUTE - 1) * 60 + 30

nMaxSolutions = 5

global stats = Vector{Dict}()
global manager = LSManager(nMaxSolutions, HARD_TIMEOUT)

function MILP(m::MATH0462_project.Model, gurobi_env::Gurobi.Env, timeLimit::Union{Int,Nothing}=nothing;
    with_callback::Bool=true, with_start_solution=false, nurse_constraints=true,
    S1::Bool=true, S2::Bool=true, S3::Bool=true, S4::Bool=true, S5::Bool=true, S6::Bool=true)

    # Create a new JuMP model with Gurobi as the solver
    jm = JuMP.Model(() -> Gurobi.Optimizer(gurobi_env))
    set_optimizer_attribute(jm, "OutputFlag", 1)
    if !isnothing(timeLimit)
        set_optimizer_attribute(jm, "TimeLimit", timeLimit)
    end

    # variables hard
    @variable(jm, p_room[1:m.P, 1:m.R], Bin) # = 1 if patient p is in room r
    @variable(jm, p_present[1:m.P, 1:m.R, 1:m.D], Bin) # = 1 if patient p is present in room r on day d
    @variable(jm, r_gender[1:m.R, 1:m.D], Bin) # = 1 if room r is associated with gender 2 ("B") on day d
    @variable(jm, p_admission[1:m.P, 1:m.D], Bin) # = 1 if patient p admitted on day d
    if nurse_constraints
        @variable(jm, n_schedule[1:m.N, 1:m.R, 1:m.D], Bin) # = 1 if nurse n is scheduled in room r on day d
    end

    # variables soft
    if S1
        @variable(jm, r_max_age[1:m.R, 1:m.D] >= 0, Int) # maximum age group in room r on day d
        @variable(jm, r_min_age[1:m.R, 1:m.D] >= 0, Int) # minimum age group in room r on day d
        @variable(jm, r_delta_age[1:m.R, 1:m.D] >= 0, Int) # maximum age group difference in room r on day d
    end
    if S2 || S4
        # trouver un meilleur nom
        @variable(jm, indicator[p in 1:m.P, d in 1:m.D, 1:min(m.D - d + 1, m.patient_length_of_stay[p])], Bin) # = 1 if patient p is on its t^th day since admission on day d
    end
    if S2 || S3 || S4
        @variable(jm, patient_and_nurse[1:m.N, 1:m.P+m.O, 1:m.R, 1:m.D], Bin) # = 1 if nurse n and patient p are in room r on day d
    end
    if S2
        @variable(jm, p_skill[1:m.P, 1:m.D] >= 0, Int) # skill of patient p on day d
        @variable(jm, penalty_skill[1:m.P+m.O, 1:m.D] >= 0, Int) # the penalty due to nurse skill assigned to patient p on day d
    end
    if S3
        @variable(jm, in_charge[1:m.N, 1:m.P+m.O], Bin) # = 1 if nurse n is in charge of patient p
    end
    if S4
        @variable(jm, p_workload[1:m.P, 1:m.D] >= 0, Int) # workload of patient p on day d
        @variable(jm, p_and_n_workload[1:m.N, 1:m.P+m.O, 1:m.R, 1:m.D] >= 0, Int) # workload of nurse n due to patient p in room r on day d
        @variable(jm, n_excess_workload[1:m.N, 1:m.D] >= 0, Int) # excess workload of nurse n on day d
    end
    if S5
        @variable(jm, p_delay[1:m.P] >= 0, Int) # admission delay of patient p
    end
    if S6
        @variable(jm, unscheduled[m.optional_patients], Bin) # = 1 if optional patient p is not scheduled
    end

    # constraints
    # occupant: gender
    for o in 1:m.O, d in 1:min(m.D, m.patient_length_of_stay[m.P+o])
        @constraint(jm, r_gender[m.occupant_room[o], d] == (m.patient_gender[m.P+o] == 2))
    end

    # patient cannot be transferred from one room to another
    for p in 1:m.P
        if !(m.patient_mandatory[p]) # since H4
            @constraint(jm, sum(p_room[p, r] for r in 1:m.R) <= 1)
        end
        for d in 1:m.D
            @constraint(jm, sum(p_present[p, r, d] for r in 1:m.R) <= 1)
            @constraint(jm, sum(p_present[p, r, d] for r in 1:m.R) <= sum(p_room[p, r] for r in 1:m.R))
            for r in 1:m.R
                @constraint(jm, p_present[p, r, d] <= p_room[p, r])
            end
        end
    end

    # link admission to present
    for p in 1:m.P
        #@constraint(jm, sum(p_admission[p, d] for d in 1:m.D) <= 1) # useless since H6
        @constraint(jm, p_admission[p, 1] == sum(p_present[p, r, 1] for r in 1:m.R))
        for d in 2:m.D
            @constraint(jm, p_admission[p, d] <= sum(p_present[p, r, d] for r in 1:m.R))
        end
    end

    # continuity of presence
    for p in 1:m.P, d in 1:m.D
        d_end = min(d + m.patient_length_of_stay[p] - 1, m.D)
        @constraint(jm, (d_end - d + 1) * p_admission[p, d] <= sum(p_present[p, r, i] for r in 1:m.R, i in d:d_end))
    end

    if nurse_constraints
        # A nurse can work only on days that are in her working_days
        for n in 1:m.N, d in 1:m.D
            if !(n in m.nurse_available_on_day[d])
                for r in 1:m.R
                    @constraint(jm, n_schedule[n, r, d] == 0)
                end
            end
        end

        # single nurse to each room, for each day within the scheduling period
        for r in 1:m.R, d in 1:m.D
            @constraint(jm, sum(n_schedule[n, r, d] for n in 1:m.N) == 1)
        end
    end

    # H1: No gender mix
    for r in 1:m.R, d in 1:m.D
        @constraint(jm, r_gender[r, d] <= 1 - (sum(p_present[p, r, d] for p in 1:m.P if m.patient_gender[p] == 1) + sum(m.o_present[o, r, d] for o in 1:m.O if m.patient_gender[m.P+o] == 1)) / m.room_capacity[r])
        @constraint(jm, r_gender[r, d] >= (sum(p_present[p, r, d] for p in 1:m.P if m.patient_gender[p] == 2) + sum(m.o_present[o, r, d] for o in 1:m.O if m.patient_gender[m.P+o] == 2)) / m.room_capacity[r])
    end

    # H2: Compatible rooms
    for p in 1:m.P, r in m.patient_invalid_rooms[p]
        @constraint(jm, p_room[p, r] == 0)
    end

    # H3: Room capacity
    for r in 1:m.R, d in 1:m.D
        @constraint(jm, sum(p_present[p, r, d] for p in 1:m.P) + sum(m.o_present[o, r, d] for o in 1:m.O) <= m.room_capacity[r])
    end

    # H4: Mandatory versus optional patients
    for p in m.mandatory_patients
        @constraint(jm, sum(p_room[p, r] for r in 1:m.R) == 1)
    end

    # H5: Admission day
    for p in m.mandatory_patients
        @constraint(jm, sum(p_admission[p, d] for d in m.patient_release_day[p]:m.patient_due_day[p]) == 1)
        @constraint(jm, sum(p_admission[p, d] for d in 1:m.patient_release_day[p]-1) == 0)
        @constraint(jm, sum(p_admission[p, d] for d in m.patient_due_day[p]+1:m.D) == 0)
    end
    for p in m.optional_patients
        @constraint(jm, sum(p_admission[p, d] for d in m.patient_release_day[p]:m.D) <= 1)
        @constraint(jm, sum(p_admission[p, d] for d in 1:m.patient_release_day[p]-1) == 0)
    end

    # S1: Age groups
    if S1
        for r in 1:m.R, d in 1:m.D
            for p in 1:m.P
                @constraint(jm, r_min_age[r, d] <= m.patient_age_group[p] + m.A * (1 - p_present[p, r, d]))
                @constraint(jm, r_max_age[r, d] >= m.patient_age_group[p] * p_present[p, r, d])
            end
            for o in 1:m.O
                @constraint(jm, r_min_age[r, d] <= m.patient_age_group[m.P+o] + m.A * (1 - m.o_present[o, r, d]))
                @constraint(jm, r_max_age[r, d] >= m.patient_age_group[m.P+o] * m.o_present[o, r, d])
            end
            @constraint(jm, r_delta_age[r, d] >= r_max_age[r, d] - r_min_age[r, d])
        end
        @expression(jm, s1, m.weight_room_mixed_age * sum(r_delta_age[r, d] for r in 1:m.R, d in 1:m.D))
    else
        @expression(jm, s1, 0)
    end

    if S2 || S4
        for p in 1:m.P, d in 1:m.D
            @constraint(jm, indicator[p, d, 1] == p_admission[p, d])
            for t in 2:min(m.D - d + 1, m.patient_length_of_stay[p])
                @constraint(jm, indicator[p, d, t] == indicator[p, d, t-1])
            end
        end
    end

    if (S2 || S3 || S4) && nurse_constraints
        # if a nurse n is assigned to a room r during day d, then that nurse must be in charge of the patient present in the room.
        for n in 1:m.N, p in 1:m.P, r in 1:m.R, d in 1:m.D
            # linearization constraints of n_schedule[n, r, d] * p_present[p, r, d]
            @constraint(jm, patient_and_nurse[n, p, r, d] <= n_schedule[n, r, d])
            @constraint(jm, patient_and_nurse[n, p, r, d] <= p_present[p, r, d])
            @constraint(jm, patient_and_nurse[n, p, r, d] >= n_schedule[n, r, d] + p_present[p, r, d] - 1)
        end
        for n in 1:m.N, o in 1:m.O, r in 1:m.R, d in 1:m.D
            # linearization constraints of n_schedule[n, r, d] * m.o_present[o, r, d]
            @constraint(jm, patient_and_nurse[n, m.P+o, r, d] <= n_schedule[n, r, d])
            @constraint(jm, patient_and_nurse[n, m.P+o, r, d] <= m.o_present[o, r, d])
            @constraint(jm, patient_and_nurse[n, m.P+o, r, d] >= n_schedule[n, r, d] + m.o_present[o, r, d] - 1)
        end
    end

    # S2: Minimum skill level
    if S2
        for p in 1:m.P, d in 1:m.D
            @constraint(jm, p_skill[p, d] == sum(indicator[p, d-t+1, t] * m.patient_skill[p][t] for t in 1:min(d, m.patient_length_of_stay[p])))
            @constraint(jm, penalty_skill[p, d] >= p_skill[p, d] - sum(patient_and_nurse[n, p, r, d] * m.nurse_skill_level[n] for n in 1:m.N, r in 1:m.R))
        end
        for o in 1:m.O
            last_day = min(m.D, m.patient_length_of_stay[m.P+o])
            for d in 1:last_day
                @constraint(jm, penalty_skill[m.P+o, d] >= m.patient_skill[m.P+o][d] - sum(patient_and_nurse[n, m.P+o, r, d] * m.nurse_skill_level[n] for n in 1:m.N, r in 1:m.R))
            end
            for d in last_day+1:m.D
                @constraint(jm, penalty_skill[m.P+o, d] == 0)
            end
        end
        @expression(jm, s2, m.weight_room_nurse_skill * sum(penalty_skill[p, d] for p in 1:m.P+m.O, d in 1:m.D))
    else
        @expression(jm, s2, 0)
    end

    # S3: Continuity of care
    if S3
        #= # Non-linear formulation:
        for n in 1:m.N, p in 1:m.P+m.O
            # sum / length of stay ou faire chacun individuellement ?
            @constraint(jm, in_charge[n, p] >= sum(n_schedule[n, r, d] * p_present[p, r, d] for r in 1:m.R, d in 1:m.D) / m.patient_length_of_stay[p])
        end =#
        # sum / length of stay ou faire chacun individuellement ?
        for n in 1:m.N, p in 1:m.P+m.O
            @constraint(jm, in_charge[n, p] >= sum(patient_and_nurse[n, p, r, d] for r in 1:m.R, d in 1:m.D) / m.patient_length_of_stay[p])
        end
        #= for n in 1:m.N, p in 1:m.P+m.O, r in 1:m.R, d in 1:m.D
            @constraint(jm, in_charge[n, p] >= patient_and_nurse[n, p, r, d])
        end =#
        @expression(jm, s3, m.weight_continuity_of_care * sum(in_charge[n, p] for n in 1:m.N, p in 1:m.P+m.O))
    else
        @expression(jm, s3, 0)
    end

    # S4: Maximum workload
    if S4
        for p in 1:m.P, d in 1:m.D
            @constraint(jm, p_workload[p, d] == sum(indicator[p, d-t+1, t] * m.patient_workload[p][t] for t in 1:min(d, m.patient_length_of_stay[p])))
        end
        #= # Non-linear formulation:
        for n in 1:m.N, d in 1:m.D
            @constraint(jm, n_excess_workload[n, d] >= sum(patient_and_nurse[n, p, r, d] * p_workload[p, d] for p in 1:m.P+m.O, r in 1:m.R) - m.nurse_max_load_on_day[n, d])
        end =#
        max_workload = maximum(vcat(m.patient_workload...))
        for n in 1:m.N, d in 1:m.D
            @constraint(jm, n_excess_workload[n, d] >= sum(p_and_n_workload[n, p, r, d] for p in 1:m.P+m.O, r in 1:m.R) - m.nurse_max_load_on_day[n, d])
            for p in 1:m.P, r in 1:m.R
                @constraint(jm, p_and_n_workload[n, p, r, d] <= max_workload * patient_and_nurse[n, p, r, d])
                @constraint(jm, p_and_n_workload[n, p, r, d] <= p_workload[p, d])
                @constraint(jm, p_and_n_workload[n, p, r, d] >= max_workload * patient_and_nurse[n, p, r, d] + p_workload[p, d] - max_workload)
            end
        end
        for n in 1:m.N, o in 1:m.O, r in 1:m.R
            last_day = min(m.D, m.patient_length_of_stay[m.P+o])
            for d in 1:last_day
                @constraint(jm, p_and_n_workload[n, m.P+o, r, d] <= max_workload * patient_and_nurse[n, m.P+o, r, d])
                @constraint(jm, p_and_n_workload[n, m.P+o, r, d] <= m.patient_workload[m.P+o][d])
                @constraint(jm, p_and_n_workload[n, m.P+o, r, d] >= max_workload * patient_and_nurse[n, m.P+o, r, d] + m.patient_workload[m.P+o][d] - max_workload)
            end
            for d in last_day+1:m.D
                @constraint(jm, p_and_n_workload[n, m.P+o, r, d] == 0)
            end
        end
        @expression(jm, s4, m.weight_nurse_eccessive_workload * sum(n_excess_workload[n, d] for n in 1:m.N, d in 1:m.D))
    else
        @expression(jm, s4, 0)
    end

    # S5: Admission delay
    if S5
        for p in 1:m.P
            @constraint(jm, p_delay[p] >= sum(p_admission[p, d] * d for d in 1:m.D) - m.patient_release_day[p])
        end
        @expression(jm, s5, m.weight_patient_delay * sum(p_delay[p] for p in 1:m.P))
    else
        @expression(jm, s5, 0)
    end

    # S6: Unscheduled patients
    if S6
        for p in m.optional_patients
            @constraint(jm, unscheduled[p] == 1 - sum(p_admission[p, d] for d in 1:m.D))
        end
        @expression(jm, s6, m.weight_unscheduled_optional * sum(unscheduled[p] for p in m.optional_patients; init=0))
    else
        @expression(jm, s6, 0)
    end

    if with_start_solution
        for p in 1:m.P
            for r in 1:m.R
                if m.v_patient_room[p].value == r
                    set_start_value(p_room[p, r], 1)
                else
                    set_start_value(p_room[p, r], 0)
                end
            end
            for d in 1:m.D
                if m.v_patient_admission_day[p].value == d
                    set_start_value(p_admission[p, d], 1)
                else
                    set_start_value(p_admission[p, d], 0)
                end
            end
            if S5
                set_start_value(p_delay[p], m.s_admission_delay[p].value)
            end
        end
        for r in 1:m.R, d in 1:m.D
            patient_present = m.c_patients_present_in_room[r, d].value
            for p in patient_present
                if p <= m.P
                    set_start_value(p_present[p, r, d], 1)
                end
            end
            for p in setdiff(1:m.P, patient_present)
                set_start_value(p_present[p, r, d], 0)
            end
            set_start_value(r_gender[r, d], (m.c_room_gender_n[r, d, 2].value > 0) ? 1 : 0)
        end
        if S6
            for p in m.optional_patients
                set_start_value(unscheduled[p], m.s_unscheduled[p].value)
            end
        end
    end

    # objective
    @objective(jm, Min, s1 + s2 + s3 + s4 + s5 + s6)

    function updateModel!(m, room_value, admission_value, nurse_value=[])
        for p in 1:m.P
            admission_found = false
            for d in 1:m.D
                if admission_value[p, d] > 0.5
                    update!(m.v_patient_admission_day[p], d)
                    admission_found = true
                    break
                end
            end
            if admission_found
                for r in 1:m.R
                    if room_value[p, r] > 0.5
                        update!(m.v_patient_room[p], r)
                        break
                    end
                end
            else
                update!(m.v_patient_admission_day[p], m.D + 1) # Not admitted
            end
        end
        if nurse_value != []
            for n in 1:m.N, r in 1:m.R, d in 1:m.D
                if nurse_value[n, r, d] > 0.5
                    update!(m.v_room_nurse[r, d], n)
                end
            end
        end
    end

    function callback_function(cb_data, cb_where::Cint)
        if cb_where == Gurobi.GRB_CB_MIPSOL
            m_copy = MATH0462_project.Model(m)
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            room_value = callback_value.(Ref(cb_data), p_room)
            admission_value = callback_value.(Ref(cb_data), p_admission)
            updateModel!(m_copy, room_value, admission_value)

            addWaitingFeasibleSolution!(manager, m_copy)
        end
    end

    if with_callback
        set_attribute(jm, Gurobi.CallbackFunction(), callback_function)
    end

    optimize!(jm)

    if !with_callback
        if nurse_constraints
            updateModel!(m, value.(p_room), value.(p_admission), value.(n_schedule))
        else
            updateModel!(m, value.(p_room), value.(p_admission))
        end
        return m
    end
end

function heuristics_local_search!(m::MATH0462_project.Model, gurobi_env::Gurobi.Env, ordering)
    println("heuristics_local_search")

    run_until_convergence(manager.startTime, manager.limitTime, [
        () -> LNSNurses.relax_nurses(m, gurobi_env=gurobi_env),
            () -> PatientOPT.patient_2_opt(manager.startTime, manager.limitTime, m, ordering, until_convergence=true),
            () -> PatientOPT.patient_3_opt(manager.startTime, manager.limitTime, m, ordering, until_convergence=true),
        ], after_each=() -> begin
            addSolution!(manager, m)
            yield()
        end)
    addSolution!(manager, m)
end

instance_directory = "competition_instances"
instance_prefix = "i"
solution_folder = "solutions"
function run_MILP(instances)
    gurobi_env = Gurobi.Env()
    project_root = @__DIR__
    println(project_root)
    for instance in instances
        println("MILP for instance $instance")
        start_time = time()
        instance_name = @sprintf("%s%02d.json", instance_prefix, instance)
        instance_path = joinpath(project_root, instance_directory, instance_name)
        solution_name = @sprintf("sol_%s", instance_name)

        solution_directory = joinpath(project_root, solution_folder, "MILP") # "hard_only")
        mkpath(solution_directory)
        solution_path = joinpath(solution_directory, solution_name)

        solution = MILP(MATH0462_project.Model(instance_path), gurobi_env, nothing, with_callback=false) # S1=false, S2=false, S3=false, S4=false, S5=false, S6=false)
        writeSolution!(solution, solution_path)
        println("Execution Time: $(time()-start_time)\n")
        printStats(solution)
        get_violations_scores(instance_path, solution_path)
    end
end

function run_heuristics(instances)
    gurobi_env = Gurobi.Env()
    project_root = @__DIR__
    for instance in instances
        println("Heuristics for instance $instance")
        global stats = Vector{Dict}()
        global manager = LSManager(nMaxSolutions, HARD_TIMEOUT)

        instance_name = @sprintf("%s%02d.json", instance_prefix, instance)
        instance_path = joinpath(project_root, instance_directory, instance_name)
        solution_name = @sprintf("sol_%s", instance_name)

        solution_directory = joinpath(project_root, solution_folder, "heuristic")
        mkpath(solution_directory)
        solution_path = joinpath(solution_directory, solution_name)

        m = MATH0462_project.Model(instance_path)
        m_first_available = Heuristic_first_available.heuristic_first_available(MATH0462_project.Model(m))
        if isnothing(m_first_available)
            # Find a feasible solution using only hard constraints
            m_first_available = MILP(m, gurobi_env, floor(Int, SOFT_TIMEOUT - (time() - manager.startTime)), with_callback=false, nurse_constraints=false, S1=false, S2=false, S3=false, S4=false, S5=false, S6=false)
        end
        # Find other feasible solutions using the MILP with only S5, S6 and using m_first_available as start values
        MILP(m_first_available, gurobi_env, floor(Int, SOFT_TIMEOUT - (time() - manager.startTime)), with_callback=true, with_start_solution=true, nurse_constraints=false, S1=false, S2=false, S3=false, S4=false, S5=true, S6=true)
        # Perform local search for the time left
        while time() - manager.startTime < HARD_TIMEOUT
            m = getWaitingFeasibleSolution!(manager)
            if m !== nothing
                heuristics_local_search!(MATH0462_project.Model(m), gurobi_env, randomOrdering(m))
            else
                #No more feasible solution
                break
            end
        end
        try
            solution = getBestSolution(manager)
            writeSolution!(solution, solution_path)
            printStats(solution)
            get_violations_scores(instance_path, solution_path)
        catch e
            println("No solution were found in $MINUTE minutes")
        end
    end
end
