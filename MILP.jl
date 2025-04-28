#= import Pkg
Pkg.add("JuMP")
Pkg.add("Gurobi")
Pkg.add("JSON3")
Pkg.add("Printf") =#

using JuMP, Gurobi, JSON3, Printf

include("validator.jl")

function MILP(instance_path::String, solution_path::String, S1::Bool=true, S2::Bool=true, S3::Bool=true, S4::Bool=true, S5::Bool=true, S6::Bool=true)
    instance = JSON3.read(instance_path)

    D = instance["days"]
    L = instance["skill_levels"]
    P = length(instance["patients"])
    R = length(instance["rooms"])
    A = length(instance["age_groups"])
    N = length(instance["nurses"])
    O = length(instance["occupants"])

    weight_room_mixed_age = instance["weights"]["room_mixed_age"]
    weight_room_nurse_skill = instance["weights"]["room_nurse_skill"]
    weight_continuity_of_care = instance["weights"]["continuity_of_care"]
    weight_nurse_eccessive_workload = instance["weights"]["nurse_eccessive_workload"]
    weight_patient_delay = instance["weights"]["patient_delay"]
    weight_unscheduled_optional = instance["weights"]["unscheduled_optional"]

    agegroup_id = Dict(j => i for (i, j) in enumerate(instance["age_groups"]))
    room_id = Dict(j["id"] => i for (i, j) in enumerate(instance["rooms"]))
    id_to_patient = [p["id"] for p in instance["patients"]]
    id_to_room = [r["id"] for r in instance["rooms"]]
    id_to_nurse = [n["id"] for n in instance["nurses"]]

    patients_then_occupants = [instance["patients"]; instance["occupants"]]
    optional_patients = [i for (i, patient) in enumerate(instance["patients"]) if patient["mandatory"] == false]
    mandatory_patients = [i for (i, patient) in enumerate(instance["patients"]) if patient["mandatory"] == true]

    patient_length_of_stay = [p["length_of_stay"] for p in patients_then_occupants]
    patient_workload = [p["workload_produced"] for p in patients_then_occupants]
    max_workload = maximum(vcat(patient_workload...))
    patient_skill = [p["skill_level_required"] for p in patients_then_occupants]
    patient_gender = [p["gender"] == "A" ? 1 : 2 for p in patients_then_occupants]
    patient_age_group = [agegroup_id[p["age_group"]] for p in patients_then_occupants]
    patient_invalid_rooms = [Set{Int}(room_id[r] for r in p["incompatible_room_ids"]) for p in instance["patients"]]
    patient_release_day = [p["surgery_release_day"] + 1 for p in instance["patients"]]
    patient_due_day = [get(p, "surgery_due_day", D) + 1 for p in instance["patients"]]
    patient_mandatory = [p["mandatory"] for p in instance["patients"]]

    occupant_room = [room_id[o["room_id"]] for o in instance["occupants"]]
  
    room_capacity = [r["capacity"] for r in instance["rooms"]]

    nurse_available_on_day = [Set{Int}() for i in 1:D]
    nurse_max_load_on_day = [0 for n in 1:N, i in 1:D]
    nurse_skill_level = [n["skill_level"] for n in instance["nurses"]]

    for n in 1:N, ws in instance["nurses"][n]["working_days"]
        d = ws["day"] + 1
        load = ws["max_load"]
        nurse_max_load_on_day[n, d] = load
        push!(nurse_available_on_day[d], n)
    end

    o_present = zeros(Bool, O, R, D)
    for o in 1:O, d in 1:min(D, patient_length_of_stay[P+o])
        o_present[o, occupant_room[o], d] = true
    end

    # Create a new JuMP model with Gurobi as the solver
    jm = direct_model(Gurobi.Optimizer()) #it is faster than Model(Gurobi.Optimizer) since no need to store on disk
    set_optimizer_attribute(jm, "OutputFlag", 1)
    #set_optimizer_attribute(jm, "TimeLimit", 10*60) # A terme supprimer mais ici permet de stopper tout apres 10min

    # variables hard
    @variable(jm, p_room[1:P, 1:R], Bin) # = 1 if patient p is in room r
    @variable(jm, p_present[1:P, 1:R, 1:D], Bin) # = 1 if patient p is present in room r on day d
    @variable(jm, r_gender[1:R, 1:D], Bin) # = 1 if room r is associated with gender 2 ("B") on day d
    @variable(jm, p_admission[1:P, 1:D], Bin) # = 1 if patient p admitted on day d
    @variable(jm, n_schedule[1:N, 1:R, 1:D], Bin) # = 1 if nurse n is scheduled in room r on day d
    
    # variables soft
    if S1
        @variable(jm, r_max_age[1:R, 1:D] >= 0, Int) # maximum age group in room r on day d
        @variable(jm, r_min_age[1:R, 1:D] >= 0, Int) # minimum age group in room r on day d
        @variable(jm, r_delta_age[1:R, 1:D] >= 0, Int) # maximum age group difference in room r on day d
    end
    if S2 || S4
        # trouver un meilleur nom
        @variable(jm, indicator[p in 1:P, d in 1:D, 1:min(D-d+1, patient_length_of_stay[p])], Bin) # = 1 if patient p is on its t^th day since admission on day d
    end
    if S2 || S3 || S4
        @variable(jm, patient_and_nurse[1:N, 1:P+O, 1:R, 1:D], Bin) # = 1 if nurse n and patient p are in room r on day d
    end
    if S2
        @variable(jm, p_skill[1:P, 1:D] >= 0, Int) # skill of patient p on day d
        @variable(jm, penalty_skill[1:P+O, 1:D] >= 0, Int) # the penalty due to nurse skill assigned to patient p on day d
    end
    if S3
        @variable(jm, in_charge[1:N, 1:P+O], Bin) # = 1 if nurse n is in charge of patient p
    end
    if S4
        @variable(jm, p_workload[1:P, 1:D] >= 0, Int) # workload of patient p on day d
        @variable(jm, p_and_n_workload[1:N, 1:P+O, 1:R, 1:D] >= 0, Int) # workload of nurse n due to patient p in room r on day d
        @variable(jm, n_excess_workload[1:N, 1:D] >= 0, Int) # excess workload of nurse n on day d
    end
    if S5
        @variable(jm, p_delay[1:P] >= 0) # admission delay of patient p
    end
    if S6
        @variable(jm, unscheduled[optional_patients], Bin) # = 1 if optional patient p is not scheduled
    end

    # constraints
    # occupant: gender
    for o in 1:O, d in 1:min(D, patient_length_of_stay[P+o])
        @constraint(jm, r_gender[occupant_room[o], d] == (patient_gender[P+o] == 1 ? 0 : 1))
    end

    # patient cannot be transferred from one room to another
    for p in 1:P
        if !(patient_mandatory[p]) # since H4
            @constraint(jm, sum(p_room[p, r] for r in 1:R) <= 1)
        end
        for d in 1:D
            @constraint(jm, sum(p_present[p, r, d] for r in 1:R) <= 1)
            @constraint(jm, sum(p_present[p, r, d] for r in 1:R) <= sum(p_room[p, r] for r in 1:R))
            for r in 1:R
                @constraint(jm, p_present[p, r, d] <= p_room[p, r])
            end
        end
    end

    # link admission to present
    for p in 1:P
        #@constraint(jm, sum(p_admission[p, d] for d in 1:D) <= 1) # useless since H6
        @constraint(jm, p_admission[p, 1] == sum(p_present[p, r, 1] for r in 1:R))
        for d in 2:D
            @constraint(jm, p_admission[p, d] <= sum(p_present[p, r, d] for r in 1:R))
        end
    end

    # continuity of presence
    for p in 1:P, d in 1:D
        d_end = min(d + patient_length_of_stay[p] - 1, D)
        @constraint(jm, (d_end - d + 1) * p_admission[p, d] <= sum(p_present[p, r, i] for r in 1:R, i in d:d_end))
    end

    # A nurse can work only on days that are in her working_days
    for n in 1:N, d in 1:D
        if !(n in nurse_available_on_day[d])
            for r in 1:R
                @constraint(jm, n_schedule[n, r, d] == 0)
            end
        end
    end

    # single nurse to each room, for each day within the scheduling period
    for r in 1:R, d in 1:D
        @constraint(jm, sum(n_schedule[n, r, d] for n in 1:N) == 1)
    end

    # H1: No gender mix
    for r in 1:R, d in 1:D 
        @constraint(jm, r_gender[r, d] <= 1 - (sum(p_present[p, r, d] for p in 1:P if patient_gender[p] == 1) + sum(o_present[o, r, d] for o in 1:O if patient_gender[P+o] == 1)) / room_capacity[r])
        @constraint(jm, r_gender[r, d] >= (sum(p_present[p, r, d] for p in 1:P if patient_gender[p] == 2) + sum(o_present[o, r, d] for o in 1:O if patient_gender[P+o] == 2)) / room_capacity[r])
    end

    # H2: Compatible rooms
    for p in 1:P, r in patient_invalid_rooms[p]
        @constraint(jm, p_room[p, r] == 0)
    end

    # H3: Room capacity
    for r in 1:R, d in 1:D
        @constraint(jm, sum(p_present[p, r, d] for p in 1:P) + sum(o_present[o, r, d] for o in 1:O) <= room_capacity[r])
    end

    # H4: Mandatory versus optional patients
    for p in mandatory_patients
        @constraint(jm, sum(p_room[p, r] for r in 1:R) == 1)
    end

    # H5: Admission day
    for p in mandatory_patients
        @constraint(jm, sum(p_admission[p, d] for d in patient_release_day[p]:patient_due_day[p]) == 1)
        @constraint(jm, sum(p_admission[p, d] for d in 1:patient_release_day[p]-1) == 0)
        @constraint(jm, sum(p_admission[p, d] for d in patient_due_day[p]+1:D) == 0)
    end
    for p in optional_patients
        @constraint(jm, sum(p_admission[p, d] for d in patient_release_day[p]:D) <= 1)
        @constraint(jm, sum(p_admission[p, d] for d in 1:patient_release_day[p]-1) == 0)
    end

    # S1: Age groups
    if S1
        for r in 1:R, d in 1:D
            for p in 1:P
                @constraint(jm, r_min_age[r, d] <= patient_age_group[p] + A * (1 - p_present[p, r, d]))
                @constraint(jm, r_max_age[r, d] >= patient_age_group[p] * p_present[p, r, d])
            end
            for o in 1:O
                @constraint(jm, r_min_age[r, d] <= patient_age_group[P+o] + A * (1 - o_present[o, r, d]))
                @constraint(jm, r_max_age[r, d] >= patient_age_group[P+o] * o_present[o, r, d])
            end
            @constraint(jm, r_delta_age[r, d] >= r_max_age[r, d] - r_min_age[r, d])
        end
        @expression(jm, s1, weight_room_mixed_age * sum(r_delta_age[r, d] for r in 1:R, d in 1:D))
    else
        @expression(jm, s1, 0)
    end

    if S2 || S4
        for p in 1:P, d in 1:D
            @constraint(jm, indicator[p, d, 1] == p_admission[p, d])
            for t in 2:min(D-d+1, patient_length_of_stay[p])
                @constraint(jm, indicator[p, d, t] == indicator[p, d, t-1])
            end
        end
    end

    if S2 || S3 || S4
        # if a nurse n is assigned to a room r during day d, then that nurse must be in charge of the patient present in the room.
        for n in 1:N, p in 1:P, r in 1:R, d in 1:D
            # linearization constraints of n_schedule[n, r, d] * p_present[p, r, d]
            @constraint(jm, patient_and_nurse[n, p, r, d] <= n_schedule[n, r, d])
            @constraint(jm, patient_and_nurse[n, p, r, d] <= p_present[p, r, d])
            @constraint(jm, patient_and_nurse[n, p, r, d] >= n_schedule[n, r, d] + p_present[p, r, d] - 1)
        end
        for n in 1:N, o in 1:O, r in 1:R, d in 1:D
            # linearization constraints of n_schedule[n, r, d] * o_present[o, r, d]
            @constraint(jm, patient_and_nurse[n, P+o, r, d] <= n_schedule[n, r, d])
            @constraint(jm, patient_and_nurse[n, P+o, r, d] <= o_present[o, r, d])
            @constraint(jm, patient_and_nurse[n, P+o, r, d] >= n_schedule[n, r, d] + o_present[o, r, d] - 1)
        end
    end

    # S2: Minimum skill level
    if S2
        for p in 1:P, d in 1:D
            @constraint(jm, p_skill[p, d] == sum(indicator[p, d-t+1, t] * patient_skill[p][t] for t in 1:min(d, patient_length_of_stay[p])))
            @constraint(jm, penalty_skill[p, d] >= p_skill[p, d] - sum(patient_and_nurse[n, p, r, d] * nurse_skill_level[n] for n in 1:N, r in 1:R))
        end
        for o in 1:O, d in 1:D
            @constraint(jm, penalty_skill[P+o, d] >= patient_skill[P+o][d] - sum(patient_and_nurse[n, P+o, r, d] * nurse_skill_level[n] for n in 1:N, r in 1:R))
        end
        @expression(jm, s2, weight_room_nurse_skill * sum(penalty_skill[p, d] for p in 1:P+O, d in 1:D))
    else
        @expression(jm, s2, 0)
    end

    # S3: Continuity of care
    if S3
        #= # Non-linear formulation:
        for n in 1:N, p in 1:O+P
            # sum / length of stay ou faire chacun individuellement ?
            @constraint(jm, in_charge[n, p] >= sum(n_schedule[n, r, d] * p_present[p, r, d] for r in 1:R, d in 1:D) / patient_length_of_stay[p])
        end =#
        # sum / length of stay ou faire chacun individuellement ?
        for n in 1:N, p in 1:P+O
            @constraint(jm, in_charge[n, p] >= sum(patient_and_nurse[n, p, r, d] for r in 1:R, d in 1:D) / patient_length_of_stay[p])
        end
        #= for n in 1:N, p in 1:P+O, r in 1:R, d in 1:D
            @constraint(jm, in_charge[n, p] >= patient_and_nurse[n, p, r, d])
        end =#
        @expression(jm, s3, weight_continuity_of_care * sum(in_charge[n, p] for n in 1:N, p in 1:P+O))
    else
        @expression(jm, s3, 0)
    end

    # S4: Maximum workload
    if S4
        for p in 1:P, d in 1:D
            @constraint(jm, p_workload[p, d] == sum(indicator[p, d-t+1, t] * patient_workload[p][t] for t in 1:min(d, patient_length_of_stay[p])))
        end
        #= # Non-linear formulation:
        for n in 1:N, d in 1:D
            @constraint(jm, n_excess_workload[n, d] >= sum(patient_and_nurse[n, p, r, d] * p_workload[p, d] for p in 1:P+O, r in 1:R) - nurse_max_load_on_day[n, d])
        end =#
        for n in 1:N, d in 1:D
            @constraint(jm, n_excess_workload[n, d] >= sum(p_and_n_workload[n, p, r, d] for p in 1:P+O, r in 1:R) - nurse_max_load_on_day[n, d])
            for p in 1:P, r in 1:R
                @constraint(jm, p_and_n_workload[n, p, r, d] <= max_workload * patient_and_nurse[n, p, r, d])
                @constraint(jm, p_and_n_workload[n, p, r, d] <= p_workload[p, d])
                @constraint(jm, p_and_n_workload[n, p, r, d] >= max_workload * patient_and_nurse[n, p, r, d] + p_workload[p, d] - max_workload)
            end
            for o in 1:O, r in 1:R
                @constraint(jm, p_and_n_workload[n, P+o, r, d] <= max_workload * patient_and_nurse[n, P+o, r, d])
                @constraint(jm, p_and_n_workload[n, P+o, r, d] <= patient_workload[P+o][d])
                @constraint(jm, p_and_n_workload[n, P+o, r, d] >= max_workload * patient_and_nurse[n, P+o, r, d] + patient_workload[P+o][d] - max_workload)
            end
        end
        @expression(jm, s4, weight_nurse_eccessive_workload * sum(n_excess_workload[n, d] for n in 1:N, d in 1:D))
    else
        @expression(jm, s4, 0)
    end

    # S5: Admission delay
    if S5
        for p in 1:P
            @constraint(jm, p_delay[p] >= sum(p_admission[p, d] * d for d in 1:D) - patient_release_day[p])
        end
        @expression(jm, s5, weight_patient_delay * sum(p_delay[p] for p in 1:P))
    else
        @expression(jm, s5, 0)
    end

    # S6: Unscheduled patients
    if S6
        for p in optional_patients
            @constraint(jm, unscheduled[p] == 1 - sum(p_admission[p, d] for d in 1:D))
        end
        @expression(jm, s6, weight_unscheduled_optional * sum(unscheduled[p] for p in optional_patients; init=0))
    else
        @expression(jm, s6, 0)
    end

    # objective
    @objective(jm, Min, s1 + s2 + s3 + s4 + s5 + s6)

    optimize!(jm)
    if termination_status(jm) == MOI.OPTIMAL || primal_status(jm) == FEASIBLE_POINT
        if termination_status(jm) == MOI.TIME_LIMIT
            println("Time limit reached in MILP\n")
        end

        solution = Dict("patients" => [], "nurses" => [])
        for p in 1:P
            patient_data = Dict{String, Union{Int64, String}}(
                "id" => id_to_patient[p],
                "admission_day" => "none"
            )

            admission_found = false
            for d in 1:D
                if value(p_admission[p, d]) > 0.5
                    patient_data["admission_day"] = d-1
                    admission_found = true
                    break
                end
            end

            if admission_found
                for r in 1:R
                    if value(p_room[p, r]) > 0.5
                        patient_data["room"] = id_to_room[r]
                        break
                    end
                end
            end

            push!(solution["patients"], patient_data)
        end

        for n in 1:N
            assignments = []
            for d in 1:D
                assigned_rooms = []
                for r in 1:R
                    if value(n_schedule[n, r, d]) > 0.5
                        push!(assigned_rooms, id_to_room[r])
                    end
                end
                if !isempty(assigned_rooms)
                    push!(assignments, Dict(
                        "day" => d-1,
                        "rooms" => assigned_rooms
                    ))
                end
            end
            push!(solution["nurses"], Dict(
                "id" => id_to_nurse[n],
                "assignments" => assignments
            ))
        end

        open(solution_path, "w") do io
            JSON3.pretty(io, solution)
        end

        println("Verification that what we find is like in the validator")
        e_total_val, s_total_val, s1_val, s2_val, s3_val, s4_val, s5_val, s6_val = get_violations_scores(instance_path, solution_path)

        if e_total_val != 0
            @printf "Violations : %7d\n" e_total_val
        end
        @printf "Objective: %7d vs validator: %7d same? %s \n" objective_value(jm) s_total_val objective_value(jm)==s_total_val
        @printf "\tS1 (age groups)      %7d vs validator: %7d same? %s \n" value(s1) s1_val value(s1)==s1_val
        @printf "\tS2 (skill level)     %7d vs validator: %7d same? %s \n" value(s2) s2_val value(s2)==s2_val
        @printf "\tS3 (continuity)      %7d vs validator: %7d same? %s \n" value(s3) s3_val value(s3)==s3_val
        @printf "\tS4 (nurse workload)  %7d vs validator: %7d same? %s \n" value(s4) s4_val value(s4)==s4_val
        @printf "\tS5 (admission delay) %7d vs validator: %7d same? %s \n" value(s5) s5_val value(s5)==s5_val
        @printf "\tS6 (unscheduled)     %7d vs validator: %7d same? %s \n" value(s6) s6_val value(s6)==s6_val

        open("times.txt", "a") do f
            write(f, match(r".*\\(.*\\.*)$", solution_path).captures[1] *
                    " " * join(["S$i" for (i, val) in enumerate([S1, S2, S3, S4, S5, S6]) if val], ",") *
                    " objective $(objective_value(jm))" *
                    " validator obj $s_total_val" *
                    " optimal ? " * string(termination_status(jm) == MOI.OPTIMAL) *
                    " Gap $(relative_gap(jm))" *
                    " Solver Time $(solve_time(jm))")
        end

        return true
    else
        println("No feasible solution found in MILP")
        open("times.txt", "a") do f
            write(f, match(r".*\\(.*\\.*)$", solution_path).captures[1] *
                    " " * join(["S$i" for (i, val) in enumerate([S1, S2, S3, S4, S5, S6]) if val], ",") *
                    " No solution" *
                    " Solver Time $(solve_time(jm))")
        end
        return false
    end
end

test_folder = false
instance_directory, instance_prefix = test_folder ? ("test_instances", "test") : ("competition_instances", "i")
solution_folder = "solutions"
for i in 1:20
    instance_name = @sprintf("%s%02d.json", instance_prefix, i)
    instance_path = joinpath(@__DIR__, instance_directory, instance_name)
    solution_name = @sprintf("sol_%s", instance_name)

    solution_directory = joinpath(@__DIR__, solution_folder, "MILP") #, "hard_only")
    mkpath(solution_directory)
    solution_path = joinpath(solution_directory, solution_name)

    start_time = time()
    MILP(instance_path, solution_path) #, false, false, false, false, false, false)
    open("times.txt", "a") do f
        write(f, " - Time: $(time()-start_time)\n")
    end

    #= solution_directory = joinpath(@__DIR__, solution_folder, "heuristic_1")
    mkpath(solution_directory)
    solution_path = joinpath(solution_directory, solution_name)
    heuristic_1(instance_path, solution_path)
    get_violations_scores(instance_path, solution_path) =#
end
