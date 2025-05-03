using Printf, JSON3

function get_violations_scores(instance_path::String, solution_path::String)
    instance = JSON3.read(instance_path)
    solution = JSON3.read(solution_path)

    D = instance["days"]
    L = instance["skill_levels"]
    P = length(instance["patients"])
    R = length(instance["rooms"])
    A = length(instance["age_groups"])
    N = length(instance["nurses"])
    O = length(instance["occupants"])

    agegroup_id = Dict(j => i for (i, j) in enumerate(instance["age_groups"]))
    room_id = Dict(j["id"] => i for (i, j) in enumerate(instance["rooms"]))

    patients_then_occupants = [instance["patients"]; instance["occupants"]]
    patient_length_of_stay = [p["length_of_stay"] for p in patients_then_occupants]
    patient_workload = [p["workload_produced"] for p in patients_then_occupants]
    patient_skill = [p["skill_level_required"] for p in patients_then_occupants]
    patient_gender = [p["gender"] == "A" ? 1 : 2 for p in patients_then_occupants]
    patient_age_group = [agegroup_id[p["age_group"]] for p in patients_then_occupants]
    patient_invalid_rooms = [Set{Int}(room_id[r] for r in p["incompatible_room_ids"]) for p in instance["patients"]]
    patient_release_day = [p["surgery_release_day"] + 1 for p in instance["patients"]]
    patient_due_day = [get(p, "surgery_due_day", D) + 1 for p in instance["patients"]]
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

    patient_base_room = [1 for p in 1:P]
    for p in 1:P
        while patient_base_room[p] in patient_invalid_rooms[p]
            patient_base_room[p] += 1
        end
    end

    # DECISIONS VARIABLES
    # Arrival date of each patient and occupants. D+1 means no admission.
    v_patient_arrival_day = [p["admission_day"] != "none" ? p["admission_day"]+1 : D+1 for p in solution["patients"]]
    v_patient_arrival_day = [v_patient_arrival_day; [1 for o in instance["occupants"]]]
    # Room of each patient and occupants.
    v_patient_room = [haskey(p, "room") ? parse(Int, p["room"][2:end])+1 : patient_base_room[parse(Int, p["id"][2:end])+1] for p in solution["patients"]]
    v_patient_room = [v_patient_room; [parse(Int, o["room_id"][2:end]) for o in instance["occupants"]]]
    # Nurse of each room
    v_room_nurse = fill(0, R, D)
    for n in solution["nurses"]
        nurse_idx = parse(Int, n["id"][2:end])
        for a in n["assignments"], r in a["rooms"]
            v_room_nurse[parse(Int, r[2:end])+1, a["day"]+1] = nurse_idx + 1
        end
    end

    # COMPUTED VARIABLES: Compute various information about the solution
    # Count the workload of a room for each room, day
    c_room_workload = [0 for i in 1:R, j in 1:D]
    # Count the workload of a nurse for each day
    c_nurse_workload = [0 for i in 1:N, j in 1:D]
    # Count the number of patients needing a given skill for each room, day, skill
    c_room_skill = [0 for i in 1:R, j in 1:D,  k in 1:L]
    # Counts the number of patients of each gender in each room for each day
    c_room_gender_n = [0 for i in 1:R, j in 1:D, g in 1:2]
    # Counts the number of patients of each age group in each room for each day
    c_room_age_n = [0 for i in 1:R, j in 1:D, a in 1:A]
    # Set of patients present in each room for each day
    c_patients_present_in_room = [Set() for i in 1:R, j in 1:D]
    # Counts the number of patients in each room for each day
    c_patients_present_in_room_n = [0 for i in 1:R, j in 1:D]
    # Nurses per patient
    c_patient_nurses = [Set() for i in 1:P+O]

    for p in solution["patients"]
        if p["admission_day"] == "none"
            continue
        end
        patient_index = parse(Int, p["id"][2:end]) + 1
        admission_day = p["admission_day"] + 1        
        room_index = parse(Int, p["room"][2:end]) + 1
        for day in 1:patient_length_of_stay[patient_index]
            if admission_day+day-1 > D
                break
            end
            nurse = v_room_nurse[room_index, admission_day+day-1]
            c_room_workload[room_index, admission_day+day-1] += patient_workload[patient_index][day]
            c_nurse_workload[nurse, admission_day+day-1] += patient_workload[patient_index][day]
            if patient_skill[patient_index][day] != 0
                c_room_skill[room_index, admission_day+day-1, patient_skill[patient_index][day]] += 1
            end
            c_room_gender_n[room_index, admission_day+day-1, patient_gender[patient_index]] += 1
            c_room_age_n[room_index, admission_day+day-1, patient_age_group[patient_index]] += 1
            push!(c_patients_present_in_room[room_index, admission_day+day-1], patient_index)
            c_patients_present_in_room_n[room_index, admission_day+day-1] += 1
            push!(c_patient_nurses[patient_index], nurse)
        end
    end
    for o in 1:O
        room_index = occupant_room[o]
        for day in 1:patient_length_of_stay[P+o]
            nurse = v_room_nurse[room_index, day]
            c_room_workload[room_index, day] += patient_workload[P+o][day]
            c_nurse_workload[nurse, day] += patient_workload[P+o][day]
            if patient_skill[P+o][day] != 0
                c_room_skill[room_index, day, patient_skill[P+o][day]] += 1
            end
            c_room_age_n[room_index, day, patient_age_group[P+o]] += 1
            push!(c_patients_present_in_room[room_index, day], P+o)
            c_patients_present_in_room_n[room_index, day] += 1
            push!(c_patient_nurses[P+o], nurse)
            c_room_gender_n[room_index, day, instance["occupants"][o]["gender"] == "A" ? 1 : 2] += 1
        end
    end
    
    # HARD CONSTRAINTS
    # H1
    e_room_gender_mix = [
        (c_room_gender_n[r, d, 1] > 0 && c_room_gender_n[r, d, 2] > 0) ? 1 : 0
        for r in 1:R, d in 1:D
    ]
    # H2
    e_patient_invalid_room = [v_patient_room[p] in patient_invalid_rooms[p] ? 1 : 0 for p in 1:P]
    e_patient_invalid_room_2 = [v_patient_room[p] <= 0 || v_patient_room[p] > R ? 1 : 0 for p in 1:P]
    # H3
    e_room_excess_capacity = [max(c_patients_present_in_room_n[r, d] - room_capacity[r], 0) for r in 1:R, d in 1:D]
    # H4
    e_patient_mandatory = [v_patient_arrival_day[p] == D + 1 ? 1 : 0 for p in 1:P if get(instance["patients"][p], "mandatory", false)]
    # H5
    e_patient_admission_day_incorrect = [patient_release_day[p] <= v_patient_arrival_day[p] <= patient_due_day[p] ? 0 : 1 for p in 1:P]
    # H6
    e_room_invalid_nurse = [in(v_room_nurse[r, d], nurse_available_on_day[d]) ? 0 : 1 for r in 1:R, d in 1:D]

    # Sum of violations
    e_patient_mandatory_n = sum(e_patient_mandatory)
    e_patient_invalid_room_n = sum(e_patient_invalid_room)
    e_patient_invalid_room_2_n = sum(e_patient_invalid_room_2)
    e_room_excess_capacity_n = sum(vec(e_room_excess_capacity))
    e_room_gender_mix_n = sum(vec(e_room_gender_mix))
    e_patient_admission_day_incorrect_n = sum(e_patient_admission_day_incorrect)
    e_room_invalid_nurse_n = sum(vec(e_room_invalid_nurse))
    e_total = sum([
        e_patient_mandatory_n * 10000000,
        e_patient_invalid_room_n,
        e_patient_invalid_room_2_n,
        e_room_excess_capacity_n,
        e_room_gender_mix_n,
        e_patient_admission_day_incorrect_n * 10000,
        e_room_invalid_nurse_n
    ])
    open(replace(solution_path, ".json" => "") * "_score.txt", "w") do f
        println(f, "Violations: ", e_total)
        if e_patient_mandatory_n != 0
            println(f, "\tMandatory patients not scheduled: ", e_patient_mandatory_n)
        end
        if e_patient_invalid_room_n != 0
            println(f, "\tPatients with invalid room: ", e_patient_invalid_room_n)
        end
        if e_room_excess_capacity_n != 0
            println(f, "\tRooms exceeding capacity: ", e_room_excess_capacity_n)
        end
        if e_room_gender_mix_n != 0
            println(f, "\tRooms with mixed gender: ", e_room_gender_mix_n)
        end
        if e_patient_admission_day_incorrect_n != 0
            println(f, "\tPatients with invalid admission day: ", e_patient_admission_day_incorrect_n)
        end
        if e_room_invalid_nurse_n != 0
            println(f, "\tRooms with invalid nurses: ", e_room_invalid_nurse_n)
        end
    end

    # SOFT CONSTRAINTS
    # S1
    s_room_delta_age = [0 for i in 1:R, j in 1:D]
    for r in 1:R, d in 1:D
        mina = A
        maxa = 1
        for a in 1:A
            if c_room_age_n[r, d, a] > 0
                mina = min(mina, a)
                maxa = max(maxa, a)
            end
        end
        s_room_delta_age[r, d] = max(maxa - mina, 0)
    end
    s_room_delta_age_total = sum(vec(s_room_delta_age))
    s1 = s_room_delta_age_total * instance["weights"]["room_mixed_age"]
    #S2
    s_room_skill_unreached = [
        sum(
            max((k - nurse_skill_level[v_room_nurse[r, d]]) * c_room_skill[r, d, k], 0) 
        )
        for r in 1:R, d in 1:D, k in 1:L
    ]


    s_room_skill_unreached_total = sum(vec(s_room_skill_unreached))
    s2 = s_room_skill_unreached_total * instance["weights"]["room_nurse_skill"]
    #S3
    s_patient_continuity = [length(c_patient_nurses[p]) for p in 1:P+O]

    s_patient_continuity_total = sum(vec(s_patient_continuity))
    s3 = s_patient_continuity_total * instance["weights"]["continuity_of_care"]
    #S4
    s_nurse_overtime = [max(c_nurse_workload[n, d] - nurse_max_load_on_day[n, d], 0) for n in 1:N, d in 1:D]
    s_nurse_overtime_total = sum(vec(s_nurse_overtime))
    s4 = s_nurse_overtime_total * instance["weights"]["nurse_eccessive_workload"]
    #S5
    s_admission_delay = [v == D + 1 ? 0 : v - patient_release_day[p] for p in 1:P for v in v_patient_arrival_day[p]]
    s_admission_delay_total = sum(vec(s_admission_delay))
    s5 = s_admission_delay_total * instance["weights"]["patient_delay"]
    #S6
    s_unscheduled = [v_patient_arrival_day[p] == D + 1 ? 1 : 0 for p in 1:P]
    s_unscheduled_total = sum(vec(s_unscheduled))
    s6 = s_unscheduled_total * instance["weights"]["unscheduled_optional"]

    s_total = sum([s1, s2, s3, s4, s5, s6])

    # Modifier ici pour avoir les scores de sauvegarder pour comparer facilement
    open(replace(solution_path, ".json" => "") * "_score.txt", "a") do f
        println(f, "Objective: ", s_total)
        @printf f "\tS1 (age groups)      %7d (%5d x %5d)\n" s1 s_room_delta_age_total instance["weights"]["room_mixed_age"]
        @printf f "\tS2 (skill level)     %7d (%5d x %5d)\n" s2 s_room_skill_unreached_total instance["weights"]["room_nurse_skill"]
        @printf f "\tS3 (continuity)      %7d (%5d x %5d)\n" s3 s_patient_continuity_total instance["weights"]["continuity_of_care"]
        @printf f "\tS4 (nurse workload)  %7d (%5d x %5d)\n" s4 s_nurse_overtime_total instance["weights"]["nurse_eccessive_workload"]
        @printf f "\tS5 (admission delay) %7d (%5d x %5d)\n" s5 s_admission_delay_total instance["weights"]["patient_delay"]
        @printf f "\tS6 (unscheduled)     %7d (%5d x %5d)\n" s6 s_unscheduled_total instance["weights"]["unscheduled_optional"]
    end

    # Modifier pour avoir le detail de chaque soft
    #return [e_total, s_total]
    return e_total, s_total, s1, s2, s3, s4, s5, s6
end