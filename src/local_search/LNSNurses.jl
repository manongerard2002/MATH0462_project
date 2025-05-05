"""
Relax everything for the Nurses
"""
module LNSNurses
using ..Reversibles
using JuMP, Gurobi
using ..LocalSearchUtils

function relax_nurses(m; gurobi_env::Gurobi.Env)
    # penalty for S2: Minimum skill level: the penalty represents the maximum between 0 and the difference between the patient and nurse skill level
    skill_level_penalty = [0 for _ in 1:m.R, _ in 1:m.D, _ in 1:m.L]
    for r in 1:m.R, d in 1:m.D, k in 1:m.L-1 #-1 since no difference of skill at max skill => no penalty
        skill_level_penalty[r, d, k] = sum(max((p_skill - k) * m.c_room_skill[r, d, p_skill].value, 0) for p_skill in 1:m.L)
    end
    nurse_skill_penalty = [0 for _ in 1:m.N, _ in 1:m.R, _ in 1:m.D]
    for n in 1:m.N, r in 1:m.R, d in 1:m.D
        nurse_skill_penalty[n, r, d] = skill_level_penalty[r, d, m.nurse_skill_level[n]]
    end
    old_objective = m.s2.value + m.s3.value + m.s4.value

    jm = JuMP.Model(() -> Gurobi.Optimizer(gurobi_env))
    set_optimizer_attribute(jm, "OutputFlag", 0)
    set_optimizer_attribute(jm, "TimeLimit", 15)

    # variable hard
    @variable(jm, n_schedule[1:m.N, 1:m.R, 1:m.D], Bin) # = 1 if nurse n is scheduled in room r on day d

    # variables soft
    @variable(jm, in_charge[1:m.P+m.O, 1:m.N], Bin) # = 1 if nurse n is in charge of patient p
    @variable(jm, n_excess_workload[1:m.N, 1:m.D] >= 0, Int) # excess workload of nurse n on day d

    # constraints
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

    # S3: Continuity of care: A nurse is in charge of a patient if they are both assigned to the same room on a day
    for r in 1:m.R, d in 1:m.D, p in m.c_patients_present_in_room[r, d].value, n in 1:m.N
        @constraint(jm, in_charge[p, n] >= n_schedule[n, r, d])
    end

    # S4: Maximum workload: The workload of a nurse on a day is the sum of the workload required in all the nurse's rooms that day
    for n in 1:m.N, d in 1:m.D
        @constraint(jm, n_excess_workload[n, d] >= sum(n_schedule[n, r, d] * m.c_room_workload[r, d].value for r in 1:m.R) - m.nurse_max_load_on_day[n, d])
    end

    # objective
    @objective(jm, Min, m.weight_room_nurse_skill * sum(nurse_skill_penalty[n, r, d] * n_schedule[n, r, d] for n in 1:m.N, r in 1:m.R, d in 1:m.D) +
                        m.weight_continuity_of_care * sum(in_charge[p, n] for p in 1:m.P+m.O, n in 1:m.N) +
                        m.weight_nurse_eccessive_workload * sum(n_excess_workload[n, d] for n in 1:m.N, d in 1:m.D))

    optimize!(jm)

    if !(termination_status(jm) == MOI.OPTIMAL || primal_status(jm) == FEASIBLE_POINT)
        # LNSNurses model not optimal/feasible (due to time limit)
        return false
    end
    if objective_value(jm) >= old_objective
        # LNSNurses solution not better (due to time limit or already the best)
        return false
    end

    orig_obj = m.s_total.value

    for n in 1:m.N, r in 1:m.R, d in 1:m.D
        if value(n_schedule[n, r, d]) > 0.5
            oldn = m.v_room_nurse[r, d].value
            if oldn != n
                update!(m.v_room_nurse[r, d], n)
            end
        end
    end

    new_obj = m.s_total.value

    if orig_obj > new_obj
        return true
    else
        return false
    end
end

end