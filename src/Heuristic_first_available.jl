"""
Try to assign the patients and nurse to the first room available
"""
module Heuristic_first_available
#using JuMP, Printf, DataStructures
using ..Reversibles

function heuristic_first_available(m)
    p_present = zeros(Bool, m.P + m.O, m.R, m.D)
    r_gender = zeros(Int, m.R, m.D)
    for o in m.P+1:m.P+m.O
        r = m.occupant_room[o-m.P]
        for d in 1:min(m.patient_length_of_stay[o], m.D)
            p_present[o, r, d] = true
            if r_gender[r, d] == 0
                r_gender[r, d] = m.patient_gender[o]
            end
        end
    end

    # Assign mandatory patient to the first room available during its stay
    for p in m.mandatory_patients
        patient_assigned = false
        for d in m.patient_release_day[p]:min(m.D, m.patient_due_day[p])
            if patient_assigned
                break
            end
            for r in m.patient_valid_rooms[p]
                stay_days = d:min(m.D, d + m.patient_length_of_stay[p] - 1)
                #check if there will be room for patient p over his length of stay and gender check 
                if all(sum(p_present[:, r, dt]) < m.room_capacity[r] &&
                    (r_gender[r, dt] == 0 || r_gender[r, dt] == m.patient_gender[p]) for dt in stay_days)
                    # Patient p is assigned to room r on day d for days in stay_days
                    for dt in stay_days
                        p_present[p, r, dt] = true
                        update!(m.v_patient_room[p], r)
                        r_gender[r, dt] = m.patient_gender[p]
                    end

                    patient_assigned = true
                    update!(m.v_patient_admission_day[p], d)
                    break
                end
            end
        end

        if !patient_assigned
            # Mandatory patient p has not been assigned within the scheduling period
            println("The heuristic first available was not able to produce a solution without violations")
            return nothing
        end
    end

    # Assign optional patient to the first room available during its stay, if possible
    for p in m.optional_patients
        patient_assigned = false
        for d in m.patient_release_day[p]:m.D
            if patient_assigned
                break
            end
            for r in m.patient_valid_rooms[p]
                stay_days = d:min(m.D, d + m.patient_length_of_stay[p] - 1)
                #check if there will be room for patient p over his length of stay and gender check 
                if all(sum(p_present[:, r, dt]) < m.room_capacity[r] &&
                   (r_gender[r, dt] == 0 || r_gender[r, dt] == m.patient_gender[p]) for dt in stay_days)
                    # Patient p is assigned to room r on day d for days in stay_days
                    for dt in stay_days
                        p_present[p, r, dt] = true
                        update!(m.v_patient_room[p], r)
                        r_gender[r, dt] = m.patient_gender[p]
                    end

                    patient_assigned = true
                    update!(m.v_patient_admission_day[p], d)
                    break
                end
            end
        end

        if !patient_assigned
            # Optional patient p has not been assigned within the scheduling period
            update!(m.v_patient_admission_day[p], m.D+1)
        end
    end

    nurses_workload = zeros(Int, m.N, m.D)
    room_workload = zeros(Int, m.R, m.D)
    room_max_skill_level = zeros(Int, m.R, m.D)
    n_schedule = zeros(Bool, m.N, m.R, m.D)
    # trying to respect skill level and nurse workload except if the nurse is alone in a day
    for d in 1:m.D
        for r in 1:m.R
            r_workload = 0
            min_skill_level = 0
            for patient in [p for p in 1:m.P+m.O if p_present[p, r, d]]
                current_day_of_staying = d - m.v_patient_admission_day[patient].value + 1
                r_workload += m.patient_workload[patient][current_day_of_staying]
                skill_level = m.patient_skill[patient][current_day_of_staying]
                if skill_level > min_skill_level
                    min_skill_level = skill_level
                end
            end
            room_workload[r, d] = r_workload
            room_max_skill_level[r, d] = min_skill_level
        end
        nurses_available = m.nurse_available_on_day[d]

        if length(nurses_available) == 1
            nurse = first(nurses_available)
            nurses_workload[nurse, d] = sum(room_workload[:, d])
            for r in 1:m.R
                # Assign nurse to room r at day d
                n_schedule[nurse, r, d] = true
                update!(m.v_room_nurse[r, d], nurse)
            end
        else
            for n in nurses_available
                personal_workload = 0
                for r in 1:m.R
                    if any(n_schedule[:, r, d]) #room already taken
                        continue
                    end
                    if personal_workload + room_workload[r, d] <= m.nurse_max_load_on_day[n, d] && room_max_skill_level[r, d] <= m.nurse_skill_level[n]
                        # Assign nurse n to room r at day d
                        personal_workload += room_workload[r, d]
                        n_schedule[n, r, d] = true
                        update!(m.v_room_nurse[r, d], n)
                    end
                end
                nurses_workload[n, d] = personal_workload
            end
        end
    end

    #if there are still rooms without nurses while respecting workload but no more skill level
    not_allocated_rooms_day = [(r, d) for r in 1:m.R, d in 1:m.D if !any(n_schedule[:, r, d])]
    if !isempty(not_allocated_rooms_day)
        days_with_remaining_rooms = unique([j[2] for j in not_allocated_rooms_day])
        for d in days_with_remaining_rooms
            rooms_list = [j[1] for j in not_allocated_rooms_day if j[2] == d]
            for n in m.nurse_available_on_day[d]
                personal_workload = nurses_workload[n, d]
                for r in rooms_list
                    if personal_workload + room_workload[r, d] <= m.nurse_max_load_on_day[n, d]
                        # Assign nurse n to room r at day d
                        personal_workload += room_workload[r, d]
                        n_schedule[n, r, d] = true
                        update!(m.v_room_nurse[r, d], n)
                    end
                end
                nurses_workload[n, d] = personal_workload
            end
        end
    end

    #if there are still rooms without nurses no attention on skill level nor workload
    not_allocated_rooms_day = [(r, d) for r in 1:m.R, d in 1:m.D if !any(n_schedule[:, r, d])]
    if !isempty(not_allocated_rooms_day)
        days_with_remaining_rooms = unique([j[2] for j in not_allocated_rooms_day])
        for d in days_with_remaining_rooms
            rooms_list = [j[1] for j in not_allocated_rooms_day if j[2] == d]
            for n in m.nurse_available_on_day[d]
                personal_workload = nurses_workload[n, d]
                for r in rooms_list
                    # Assign nurse n to room r at day d
                    personal_workload += room_workload[r, d]
                    n_schedule[n, r, d] = true
                    update!(m.v_room_nurse[r, d], n)
                end
                nurses_workload[n, d] = personal_workload
            end
        end
    end

    return m
end

export heuristic_first_available
end