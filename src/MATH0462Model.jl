module MATH0462Model
export Model, printStats, loadSolution!, writeSolution!
using JSON3, ..Reversibles, Printf
using DataStructures: OrderedSet


function reverseAction(reverse)
    if reverse
        return dec!, inc!, delete!, push!
    else
        return inc!, dec!, push!, delete!
    end
end

function viewOf(what::ReversibleInt, op::Function) # v => v in S ? 1 : 0
    v = ReversibleInt(op(what.value), what.trail)
    registerOnUpdate!(what, (oldv::Int, newv::Int) -> (update!(v, op(newv)); nothing))
    v
end

function viewOf(what::Union{ReversibleSet,ReversibleCountingSet}, op::Function)
    v = ReversibleInt(op(what.value), what.trail)
    registerOnPush!(what, (newset, addedv) -> (update!(v, op(newset)); nothing))
    registerOnDelete!(what, (newset, delv) -> (update!(v, op(newset)); nothing))
    v
end

function viewMul(what::ReversibleInt, weight)
    viewOf(what, v -> weight * v)
end

function biOperation(left::ReversibleInt, right::ReversibleInt, op::Function)
    v = ReversibleInt(op(left.value, right.value), left.trail)
    registerOnUpdate!(left, (oldv::Int, newv::Int) -> (update!(v, op(newv, right.value)); nothing))
    registerOnUpdate!(right, (oldv::Int, newv::Int) -> (update!(v, op(left.value, newv)); nothing))
    v
end

function respectsCond(what::ReversibleInt, cond::Function, valueIfOK=1, valueIfNot=0)
    viewOf(what, (v) -> cond(v) ? valueIfOK : valueIfNot)
end

function isEqualTo(what::ReversibleInt, value::Int, valueIfEqual=1, valueIfNotEqual=0)
    respectsCond(what, (nv) -> nv == value, valueIfEqual, valueIfNotEqual)
end

function isIn(what::ReversibleInt, value::OrderedSet{Int}, valueIfEqual=1, valueIfNotEqual=0)
    respectsCond(what, (nv) -> nv in value, valueIfEqual, valueIfNotEqual)
end

function sum(what, trail::Trail)
    initvalue = 0
    out = ReversibleInt(0, trail)
    for v in what
        initvalue += v.value
        registerOnUpdate!(v, (oldv::Int, newv::Int) -> (inc!(out, newv - oldv); nothing))
    end
    update!(out, initvalue)
    out
end

mutable struct Model
    instance::Any
    trail::Trail

    D::Int
    L::Int
    P::Int
    R::Int
    A::Int
    N::Int
    O::Int

    weight_room_mixed_age::Int
    weight_room_nurse_skill::Int
    weight_continuity_of_care::Int
    weight_nurse_eccessive_workload::Int
    weight_patient_delay::Int
    weight_unscheduled_optional::Int

    agegroup_id::Dict{String,Int}
    room_id::Dict{String,Int}
    nurse_id::Dict{String,Int}
    id_to_patient::Vector{String}
    id_to_room::Vector{String}
    id_to_nurse::Vector{String}

    optional_patients::Vector{Int}
    mandatory_patients::Vector{Int}
    patient_length_of_stay::Vector{Int}
    patient_workload::Vector{Vector{Int}}
    patient_skill::Vector{Vector{Int}}
    patient_gender::Vector{Int}
    patient_age_group::Vector{Int}
    patient_invalid_rooms::Vector{OrderedSet{Int}}
    patient_valid_rooms::Vector{OrderedSet{Int}}
    patient_release_day::Vector{Int}
    patient_due_day::Vector{Int}
    patient_mandatory::Vector{Bool}

    occupant_room::Vector{Int}

    room_capacity::Vector{Int}

    nurse_available_on_day::Vector{OrderedSet{Int}}
    nurse_max_load_on_day::Array{Int,2}
    nurse_skill_level::Vector{Int}

    o_present::Array{Bool, 3}

    # Admission day of each patient. D+1 means no admission
    v_patient_admission_day::Vector{ReversibleInt}
    # Room of each patient
    v_patient_room::Vector{ReversibleInt}
    # Nurse of each room, at each day
    v_room_nurse::Array{ReversibleInt,2}

    ###### COMPUTED VARIABLES
    ###### Compute various information about the solution
    # Count the workload of a room for each room, day
    c_room_workload::Array{ReversibleInt,2}
    # Count the workload of a nurse for each day
    c_nurse_workload::Array{ReversibleInt,2}
    # Count the number of patients needing a given skill for each room, day skill
    c_room_skill::Array{ReversibleInt,3}
    # Counts the number of patients of each gender in each room for each day
    c_room_gender_n::Array{ReversibleInt,3}
    # Counts the number of patients of each age group in each room for each day
    c_room_age_n::Array{ReversibleInt,3}
    # Set of patients present in each room for each day
    c_patients_present_in_room::Array{ReversibleSet{Int},2}
    # Counts the number of patients in each room for each day
    c_patients_present_in_room_n::Array{ReversibleInt,2}
    # Nurses per patient
    c_patient_nurses::Vector{ReversibleCountingSet{Int}}

    ###### SOFT CONSTRAINTS
    #S1
    s_room_delta_age::Array{ReversibleInt,2}
    s_room_delta_age_total::ReversibleInt
    s1::ReversibleInt
    #S2
    s_room_skill_unreached::Array{ReversibleInt,3}
    s_room_skill_unreached_total::ReversibleInt
    s2::ReversibleInt
    #S3
    s_patient_continuity::Array{ReversibleInt,1}
    s_patient_continuity_total::ReversibleInt
    s3::ReversibleInt
    #S4
    s_nurse_overtime::Array{ReversibleInt,2}
    s_nurse_overtime_total::ReversibleInt
    s4::ReversibleInt
    #S5
    s_admission_delay::Array{ReversibleInt,1}
    s_admission_delay_total::ReversibleInt
    s5::ReversibleInt
    #S6
    s_unscheduled::Array{ReversibleInt,1}
    s_unscheduled_total::ReversibleInt
    s6::ReversibleInt

    s_total::ReversibleInt


    ###### HARD CONSTRAINTS
    ###### as violations
    # H1
    e_room_gender_mix::Array{ReversibleInt,2}
    # H2
    e_patient_invalid_room::Array{ReversibleInt,1}
    e_patient_invalid_room_2::Array{ReversibleInt,1}
    # H3
    e_room_excess_capacity::Array{ReversibleInt,2}
    # H4
    e_patient_mandatory::Array{ReversibleInt,1}
    # H5
    e_patient_admission_day_incorrect::Array{ReversibleInt,1}
    # H6 (bonus)
    e_room_invalid_nurse::Array{ReversibleInt,2}

    # Sum of violations
    e_patient_mandatory_n::ReversibleInt
    e_patient_invalid_room_n::ReversibleInt
    e_patient_invalid_room_2_n::ReversibleInt
    e_room_excess_capacity_n::ReversibleInt
    e_room_gender_mix_n::ReversibleInt
    e_patient_admission_day_incorrect_n::ReversibleInt
    e_room_invalid_nurse_n::ReversibleInt
    e_total::ReversibleInt

    function Model(instancepath::String)
        return Model(JSON3.read(instancepath))
    end

    function Model(from::Model)
        m = Model(from.instance)

        for p in 1:m.P+m.O
            update!(m.v_patient_admission_day[p], from.v_patient_admission_day[p].value)
            update!(m.v_patient_room[p], from.v_patient_room[p].value)
        end

        for r in 1:m.R, d in 1:m.D
            update!(m.v_room_nurse[r, d], from.v_room_nurse[r, d].value)
        end

        return m
    end

    function Model(instance::JSON3.Object)
        m = new()
        m.instance = instance

        m.trail = Trail()

        m.D = instance["days"]
        m.L = instance["skill_levels"]
        m.P = length(instance["patients"])
        m.R = length(instance["rooms"])
        m.A = length(instance["age_groups"])
        m.N = length(instance["nurses"])
        m.O = length(instance["occupants"])

        m.weight_room_mixed_age = instance["weights"]["room_mixed_age"]
        m.weight_room_nurse_skill = instance["weights"]["room_nurse_skill"]
        m.weight_continuity_of_care = instance["weights"]["continuity_of_care"]
        m.weight_nurse_eccessive_workload = instance["weights"]["nurse_eccessive_workload"]
        m.weight_patient_delay = instance["weights"]["patient_delay"]
        m.weight_unscheduled_optional = instance["weights"]["unscheduled_optional"]

        m.agegroup_id = Dict(j => i for (i, j) in enumerate(instance["age_groups"]))
        m.room_id = Dict(j["id"] => i for (i, j) in enumerate(instance["rooms"]))
        m.nurse_id = Dict(j["id"] => i for (i, j) in enumerate(instance["nurses"]))
        m.id_to_patient = [p["id"] for p in instance["patients"]]
        m.id_to_room = [r["id"] for r in instance["rooms"]]
        m.id_to_nurse = [n["id"] for n in instance["nurses"]]

        patients_then_occupants = [instance["patients"]; instance["occupants"]]
        m.optional_patients = [i for (i, patient) in enumerate(instance["patients"]) if !patient["mandatory"]]
        m.mandatory_patients = [i for (i, patient) in enumerate(instance["patients"]) if patient["mandatory"]]

        m.patient_length_of_stay = [p["length_of_stay"] for p in patients_then_occupants]
        m.patient_workload = [p["workload_produced"] for p in patients_then_occupants]
        m.patient_skill = [p["skill_level_required"] .+ 1 for p in patients_then_occupants]
        m.patient_gender = [p["gender"] == "A" ? 1 : 2 for p in patients_then_occupants]
        m.patient_age_group = [m.agegroup_id[p["age_group"]] for p in patients_then_occupants]
        m.patient_invalid_rooms = [OrderedSet(m.room_id[r] for r in p["incompatible_room_ids"]) for p in instance["patients"]]
        m.patient_valid_rooms = [setdiff(OrderedSet(1:m.R), m.patient_invalid_rooms[p]) for p in 1:m.P]
        m.patient_release_day = [p["surgery_release_day"] + 1 for p in instance["patients"]]
        m.patient_due_day = [get(p, "surgery_due_day", m.D) + 1 for p in instance["patients"]]
        m.patient_mandatory = [p["mandatory"] for p in instance["patients"]]
        m.occupant_room = [m.room_id[o["room_id"]] for o in instance["occupants"]]

        m.room_capacity = [r["capacity"] for r in instance["rooms"]]

        m.nurse_available_on_day = [OrderedSet{Int}() for i in 1:m.D]
        m.nurse_max_load_on_day = [0 for n in 1:m.N, i in 1:m.D]
        m.nurse_skill_level = [n["skill_level"] .+ 1 for n in instance["nurses"]]

        for n in 1:m.N, ws in instance["nurses"][n]["working_days"]
            d = ws["day"] + 1
            load = ws["max_load"]
            m.nurse_max_load_on_day[n, d] = load
            push!(m.nurse_available_on_day[d], n)
        end

        m.o_present = zeros(Bool, m.O, m.R, m.D)
        for o in 1:m.O, d in 1:min(m.D, m.patient_length_of_stay[m.P+o])
            m.o_present[o, m.occupant_room[o], d] = true
        end

        ###### DECISIONS VARIABLES
        ###### every c_, e_, s_ variables are computed from these v_ variables automatically
        # Arrival date of each patient and occupants. D+1 means no admission
        m.v_patient_admission_day = [ReversibleInt(p <= m.P ? rand(m.patient_release_day[p]:m.patient_due_day[p]) : 1, m.trail) for p in 1:m.P+m.O]
        # Room of each patient and occupants
        m.v_patient_room = [ReversibleInt(p <= m.P ? rand(m.patient_valid_rooms[p]) : m.occupant_room[p-m.P], m.trail) for p in 1:m.P+m.O]
        # Nurse of each room
        m.v_room_nurse = [ReversibleInt(rand(m.nurse_available_on_day[d]), m.trail) for r in 1:m.R, d in 1:m.D]

        ###### COMPUTED VARIABLES
        ###### Compute various information about the solution
        # Count the workload of a room for each room, day
        m.c_room_workload = [ReversibleInt(0, m.trail) for i in 1:m.R, j in 1:m.D]
        # Count the workload of a nurse for each day
        m.c_nurse_workload = [ReversibleInt(0, m.trail) for i in 1:m.N, j in 1:m.D]
        # Count the number of patients needing a given skill for each room, day, skill
        m.c_room_skill = [ReversibleInt(0, m.trail) for i in 1:m.R, j in 1:m.D, k in 1:m.L]
        # Counts the number of patients of each gender in each room for each day
        m.c_room_gender_n = [ReversibleInt(0, m.trail) for i in 1:m.R, j in 1:m.D, g in 1:2]
        # Counts the number of patients of each age group in each room for each day
        m.c_room_age_n = [ReversibleInt(0, m.trail) for i in 1:m.R, j in 1:m.D, a in 1:m.A]
        # Set of patients present in each room for each day
        m.c_patients_present_in_room = [ReversibleSet{Int}(m.trail) for i in 1:m.R, j in 1:m.D]
        # Counts the number of patients in each room for each day
        m.c_patients_present_in_room_n = [ReversibleInt(0, m.trail) for i in 1:m.R, j in 1:m.D]
        # Nurses per patient
        m.c_patient_nurses = [ReversibleCountingSet{Int}(m.trail) for i in 1:m.P+m.O]

        ###### HARD CONSTRAINTS
        ###### as violations
        # H1
        m.e_room_gender_mix = [biOperation(m.c_room_gender_n[r, d, 1], m.c_room_gender_n[r, d, 2], (a, b) -> (a > 0 && b > 0) ? 1 : 0) for r in 1:m.R, d in 1:m.D]
        # H2
        m.e_patient_invalid_room = [isIn(m.v_patient_room[p], m.patient_invalid_rooms[p]) for p in 1:m.P]
        m.e_patient_invalid_room_2 = [respectsCond(m.v_patient_room[p], v -> v <= 0 || v > m.R) for p in 1:m.P]
        # H3
        m.e_room_excess_capacity = [viewOf(m.c_patients_present_in_room_n[r, d], v -> max(v - m.room_capacity[r], 0)) for r in 1:m.R, d in 1:m.D]
        # H4
        m.e_patient_mandatory = [isEqualTo(m.v_patient_admission_day[p], m.D + 1, 1, 0) for p in 1:m.P if m.patient_mandatory[p]]
        # H5
        m.e_patient_admission_day_incorrect = [respectsCond(m.v_patient_admission_day[p], v -> m.patient_release_day[p] <= v <= m.patient_due_day[p], 0, 1) for p in 1:m.P]
        # H6 (bonus since not present in MILP)
        m.e_room_invalid_nurse = [isIn(m.v_room_nurse[r, d], m.nurse_available_on_day[d], 0, 1) for r in 1:m.R, d in 1:m.D]

        # Sum of violations
        m.e_patient_mandatory_n = sum(m.e_patient_mandatory, m.trail)
        m.e_patient_invalid_room_n = sum(m.e_patient_invalid_room, m.trail)
        m.e_patient_invalid_room_2_n = sum(m.e_patient_invalid_room_2, m.trail)
        m.e_room_excess_capacity_n = sum(vec(m.e_room_excess_capacity), m.trail)
        m.e_room_gender_mix_n = sum(vec(m.e_room_gender_mix), m.trail)
        m.e_patient_admission_day_incorrect_n = sum(m.e_patient_admission_day_incorrect, m.trail)
        m.e_room_invalid_nurse_n = sum(vec(m.e_room_invalid_nurse), m.trail)
        m.e_total = sum([
                viewMul(m.e_patient_mandatory_n, 10000000),
                m.e_patient_invalid_room_n,
                m.e_patient_invalid_room_2_n,
                m.e_room_excess_capacity_n,
                m.e_room_gender_mix_n,
                viewMul(m.e_patient_admission_day_incorrect_n, 10000),
                m.e_room_invalid_nurse_n
            ], m.trail)

        ###### SOFT CONSTRAINTS
        #S1
        m.s_room_delta_age = [ReversibleInt(0, m.trail) for i in 1:m.R, j in 1:m.D]
        m.s_room_delta_age_total = sum(vec(m.s_room_delta_age), m.trail)
        m.s1 = viewMul(m.s_room_delta_age_total, m.weight_room_mixed_age)
        #S2
        m.s_room_skill_unreached = [biOperation(m.c_room_skill[r, d, k], m.v_room_nurse[r, d], (n_at_skill, n) -> max((k - m.nurse_skill_level[n]) * n_at_skill, 0))
                                    for r in 1:m.R, d in 1:m.D, k in 1:m.L]
        m.s_room_skill_unreached_total = sum(vec(m.s_room_skill_unreached), m.trail)
        m.s2 = viewMul(m.s_room_skill_unreached_total, m.weight_room_nurse_skill)
        #S3
        m.s_patient_continuity = [viewOf(m.c_patient_nurses[p], s -> length(s)) for p in 1:m.P+m.O]
        m.s_patient_continuity_total = sum(vec(m.s_patient_continuity), m.trail)
        m.s3 = viewMul(m.s_patient_continuity_total, m.weight_continuity_of_care)
        #S4
        m.s_nurse_overtime = [viewOf(m.c_nurse_workload[n, d], v -> max(v - m.nurse_max_load_on_day[n, d], 0)) for n in 1:m.N, d in 1:m.D]
        m.s_nurse_overtime_total = sum(vec(m.s_nurse_overtime), m.trail)
        m.s4 = viewMul(m.s_nurse_overtime_total, m.weight_nurse_eccessive_workload)
        #S5
        m.s_admission_delay = [viewOf(m.v_patient_admission_day[p], v -> v == m.D + 1 ? 0 : v - m.patient_release_day[p]) for p in 1:m.P]
        m.s_admission_delay_total = sum(vec(m.s_admission_delay), m.trail)
        m.s5 = viewMul(m.s_admission_delay_total, m.weight_patient_delay)
        #S6
        m.s_unscheduled = [isEqualTo(m.v_patient_admission_day[p], m.D + 1, 1, 0) for p in 1:m.P]
        m.s_unscheduled_total = sum(vec(m.s_unscheduled), m.trail)
        m.s6 = viewMul(m.s_unscheduled_total, m.weight_unscheduled_optional)

        m.s_total = sum([m.s1, m.s2, m.s3, m.s4, m.s5, m.s6], m.trail)

        begin
            """
            Is in charge of updating
            - c_nurse_workload
            """

            function nurseWorkloadUpdater!()
                function changed!(n, d, w, added)
                    (inc!, dec!, push!, delete!) = reverseAction(!added)
                    inc!(m.c_nurse_workload[n, d], w)
                end

                # c_room_workload[r,d] and v_room_nurse[r,d] => c_nurse_workload[n,d]
                for r in 1:m.R, d in 1:m.D
                    registerOnUpdate!(m.c_room_workload[r, d], (oldw, neww) -> (changed!(m.v_room_nurse[r, d].value, d, oldw, false); changed!(m.v_room_nurse[r, d].value, d, neww, true)))
                    registerOnUpdate!(m.v_room_nurse[r, d], (oldn, newn) -> (changed!(oldn, d, m.c_room_workload[r, d].value, false); changed!(newn, d, m.c_room_workload[r, d].value, true)))
                    changed!(m.v_room_nurse[r, d].value, d, m.c_room_workload[r, d].value, true)
                end
            end

            """
            Is in charge of updating
            - c_patient_nurses
            """

            function patientNurseUpdater!()
                function addPatient(r, d, p)
                    n = m.v_room_nurse[r, d].value
                    push!(m.c_patient_nurses[p], n)
                end

                function deletePatient(r, d, p)
                    n = m.v_room_nurse[r, d].value
                    delete!(m.c_patient_nurses[p], n)
                end

                function changeNurse(r, d, n, added)
                    (inc!, dec!, push!, delete!) = reverseAction(!added)
                    for p in m.c_patients_present_in_room[r, d]
                        push!(m.c_patient_nurses[p], n)
                    end
                end

                for r in 1:m.R, d in 1:m.D
                    registerOnPush!(m.c_patients_present_in_room[r, d], (curs, newp) -> addPatient(r, d, newp))
                    registerOnDelete!(m.c_patients_present_in_room[r, d], (curs, delp) -> deletePatient(r, d, delp))
                    registerOnUpdate!(m.v_room_nurse[r, d], (oldn, newn) -> (changeNurse(r, d, oldn, false); changeNurse(r, d, newn, true)))
                    changeNurse(r, d, m.v_room_nurse[r, d].value, true)
                end
            end

            """
            Is in charge of updating
            - c_room_workload
            - c_room_skill
            - c_patients_present_in_room
            - c_room_gender_n
            - c_room_age_n
            - c_patients_present_in_room_n
            """

            function roomWorkloadUpdater!()
                function patientChangedRW(p, day, room, added)
                    (inc!, dec!, push!, delete!) = reverseAction(!added)

                    if room == m.R + 1 || day == m.D + 1
                        return
                    end

                    for d in day:min(day + m.patient_length_of_stay[p] - 1, m.D)
                        push!(m.c_patients_present_in_room[room, d], p)
                        inc!(m.c_room_gender_n[room, d, m.patient_gender[p]], 1)
                        inc!(m.c_room_age_n[room, d, m.patient_age_group[p]], 1)
                        inc!(m.c_patients_present_in_room_n[room, d], 1)
                        inc!(m.c_room_workload[room, d], m.patient_workload[p][d-day+1])
                        inc!(m.c_room_skill[room, d, m.patient_skill[p][d-day+1]], 1)
                    end
                end

                for p in 1:m.P+m.O
                    v = m.v_patient_admission_day[p]
                    registerOnUpdate!(v, (oldv, newv) -> (patientChangedRW(p, oldv, m.v_patient_room[p].value, false); patientChangedRW(p, newv, m.v_patient_room[p].value, true)))
                    # actually init everything
                    patientChangedRW(p, m.v_patient_admission_day[p].value, m.v_patient_room[p].value, true)
                end

                for p in 1:m.P+m.O
                    v = m.v_patient_room[p]
                    registerOnUpdate!(v, (oldr, newr) -> (patientChangedRW(p, m.v_patient_admission_day[p].value, oldr, false); patientChangedRW(p, m.v_patient_admission_day[p].value, newr, true)))
                end
            end

            """
            Is in charge of updating
            - s_room_delta_age
            """

            function deltaAgeUpdater!()
                function recompute(r, d)
                    mina = m.A
                    maxa = 1
                    for a in 1:m.A
                        if m.c_room_age_n[r, d, a].value > 0
                            mina = min(mina, a)
                            maxa = max(maxa, a)
                        end
                    end
                    update!(m.s_room_delta_age[r, d], max(maxa - mina, 0))
                end

                for r in 1:m.R, d in 1:m.D
                    for a in 1:m.A
                        registerOnUpdate!(m.c_room_age_n[r, d, a], (oldn, newn) -> recompute(r, d))
                    end
                    recompute(r, d)
                end
            end
            roomWorkloadUpdater!()
            patientNurseUpdater!()
            deltaAgeUpdater!()
            nurseWorkloadUpdater!()
        end

        return m
    end
end

function printStats(m::Model)
    println("Violations: ", m.e_total.value)
    if m.e_patient_mandatory_n.value != 0
        println("\tMandatory patients not scheduled: ", m.e_patient_mandatory_n.value)
    end
    if m.e_patient_invalid_room_n.value != 0
        println("\tPatients with invalid room: ", m.e_patient_invalid_room_n.value)
    end
    if m.e_room_excess_capacity_n.value != 0
        println("\tRooms exceeding capacity: ", m.e_room_excess_capacity_n.value)
    end
    if m.e_room_gender_mix_n.value != 0
        println("\tRooms with mixed gender: ", m.e_room_gender_mix_n.value)
    end
    if m.e_patient_admission_day_incorrect_n.value != 0
        println("\tPatients with invalid admission day: ", m.e_patient_admission_day_incorrect_n.value)
    end
    if m.e_room_invalid_nurse_n.value != 0
        println("\tRooms with invalid nurses: ", m.e_room_invalid_nurse_n.value)
    end
    println("Objective: ", m.s_total.value)
    @printf "\tS1 (age groups)      %7d (%5d x %5d)\n" m.s1.value m.s_room_delta_age_total.value m.weight_room_mixed_age
    @printf "\tS2 (skill level)     %7d (%5d x %5d)\n" m.s2.value m.s_room_skill_unreached_total.value m.weight_room_nurse_skill
    @printf "\tS3 (continuity)      %7d (%5d x %5d)\n" m.s3.value m.s_patient_continuity_total.value m.weight_continuity_of_care
    @printf "\tS4 (nurse workload)  %7d (%5d x %5d)\n" m.s4.value m.s_nurse_overtime_total.value m.weight_nurse_eccessive_workload
    @printf "\tS5 (admission delay) %7d (%5d x %5d)\n" m.s5.value m.s_admission_delay_total.value m.weight_patient_delay
    @printf "\tS6 (unscheduled)     %7d (%5d x %5d)\n" m.s6.value m.s_unscheduled_total.value m.weight_unscheduled_optional
end

function writeSolution!(m::Model, path::String)
    @assert m.e_total.value == 0

    out = Dict("patients" => [], "nurses" => [])
    for p in 1:m.P
        if m.v_patient_admission_day[p].value != m.D + 1
            push!(out["patients"], Dict(
                "id" => m.id_to_patient[p],
                "admission_day" => m.v_patient_admission_day[p].value - 1,
                "room" => m.id_to_room[m.v_patient_room[p].value]
            ))
        else
            push!(out["patients"], Dict(
                "id" => m.id_to_patient[p],
                "admission_day" => "none"
            ))
        end
    end

    tmp = [Set{Int}() for n in 1:m.N, d in 1:m.D]
    for r in 1:m.R, d in 1:m.D
        n = m.v_room_nurse[r, d].value
        push!(tmp[n, d], r)
    end
    for n in 1:m.N
        assignments = []
        for d in 1:m.D
            if length(tmp[n, d]) == 0
                continue
            end
            push!(assignments, Dict(
                "rooms" => [m.id_to_room[r] for r in tmp[n, d]],
                "day" => d - 1
            ))
        end
        push!(out["nurses"], Dict(
            "id" => m.id_to_nurse[n],
            "assignments" => assignments
        ))
    end
    open(path, "w") do io
        JSON3.pretty(io, out)
    end
end
end