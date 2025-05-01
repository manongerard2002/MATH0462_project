"""import Pkg
Pkg.add("JuMP")
Pkg.add("Gurobi")
Pkg.add("JSON3")
Pkg.add("DataStructures")"""

using JuMP, Gurobi, JSON3, Printf, DataStructures

include("validator.jl")

function Heursitic(instance_path::String)
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

    room_id = Dict(j["id"] => i for (i, j) in enumerate(instance["rooms"]))

    patients_then_occupants = [instance["patients"]; instance["occupants"]]
    optional_patients = [i for (i, patient) in enumerate(instance["patients"]) if patient["mandatory"] == false]
    mandatory_patients = [i for (i, patient) in enumerate(instance["patients"]) if patient["mandatory"] == true]

    patient_length_of_stay = [p["length_of_stay"] for p in patients_then_occupants]
    patient_workload = [p["workload_produced"] for p in patients_then_occupants]
    patient_skill_level_required = [p["skill_level_required"] for p in patients_then_occupants]
    patient_gender = [p["gender"] == "A" ? 1 : 2 for p in patients_then_occupants]
    patient_invalid_rooms = [Set{Int}(room_id[r] for r in p["incompatible_room_ids"]) for p in instance["patients"]]
    patient_release_day = [p["surgery_release_day"] + 1 for p in instance["patients"]]
    patient_due_day = [get(p, "surgery_due_day", D) + 1 for p in instance["patients"]]

    occupant_room = [room_id[o["room_id"]] for o in instance["occupants"]]
  
    room_capacity = [r["capacity"] for r in instance["rooms"]]
    nurse_skill_level = [n["skill_level"] for n in instance["nurses"]]

    nurse_available_on_day_load = [
      Dict(d["day"] => d["max_load"] for d in n["working_days"])
      for n in instance["nurses"]
    ]

    p_present = zeros(Bool, P+O, R, D)
    r_gender = zeros(Int, R, D)
    for o in P+1:P+O
      r = occupant_room[o-P]
      for d in 1:min(patient_length_of_stay[o], D)
        p_present[o, r, d] = true
        if r_gender[r, d] == 0
          r_gender[r, d] = patient_gender[o]
        end
      end
    end
    
    day_entered = ones(Int, P+O)
    for p in mandatory_patients
      patient_assigned = false
      available_rooms = setdiff(OrderedSet(1:R), patient_invalid_rooms[p])
      for d in patient_release_day[p]:min(D, patient_due_day[p])
        if patient_assigned == true
          break
        end
        for r in available_rooms
          stay_days = d:min(D, d + patient_length_of_stay[p] - 1)
          #check if there will be room for patient p over his length of stay and gender check 
          if all(sum(p_present[:, r, dt]) < room_capacity[r] for dt in stay_days) &&
              all(r_gender[r, dt] == 0 || r_gender[r, dt] == patient_gender[p] for dt in stay_days)

              for dt in stay_days
                  p_present[p, r, dt] = true
                  r_gender[r, dt] = patient_gender[p]
              end

              patient_assigned = true
              day_entered[p] = d
              println("Patient $p is assigned to Room $r on Day $d for days $stay_days")
              break
          end 
        end
      end
    end

    for p in optional_patients
      patient_assigned = false
      available_rooms = setdiff(OrderedSet(1:R), patient_invalid_rooms[p])
      for d in patient_release_day[p]:D
        if patient_assigned == true
          break
        end
        for r in available_rooms
          stay_days = d:min(D, d + patient_length_of_stay[p] - 1)
          #check if there will be room for patient p over his length of stay and gender check 
          if all(sum(p_present[:, r, dt]) < room_capacity[r] for dt in stay_days) &&
              all(r_gender[r, dt] == 0 || r_gender[r, dt] == patient_gender[p] for dt in stay_days)

              for dt in stay_days
                  p_present[p, r, dt] = true
                  r_gender[r, dt] = patient_gender[p]
              end

              patient_assigned = true
              day_entered[p] = d
              println("Patient $p is assigned to Room $r on Day $d for days $stay_days")
              break
          end 
        end
      end

      if !patient_assigned
        println("Patient $p has not been assigned within the scheduling period")
      end
    end  
    nurses_workload = zeros(Int, N, D)
    room_workload = zeros(Int, R, D)
    skill_level_by_room = zeros(Int, R, D)
    assignment_matrix = zeros(Bool, N, R, D)
    # trying to respect skill level and nurse workload except if the nurse is alone in a day
    for d in 1:D
      for r in 1:R
        list_of_patients = [p for p in 1:P+O if p_present[p, r, d] == true]
        r_workload = 0
        min_skill_level = 0
        for patient in list_of_patients
          current_day_of_staying = d - day_entered[patient] + 1
          r_workload += patient_workload[patient][current_day_of_staying]
          skill_level = patient_skill_level_required[patient][current_day_of_staying]
          if skill_level > min_skill_level
            min_skill_level = skill_level
          end
        end
        room_workload[r, d] = r_workload
        skill_level_by_room[r, d] = min_skill_level
      end
      nurses_available = [n for n in 1:N if haskey(nurse_available_on_day_load[n], d-1)]
      
      if length(nurses_available) == 1
        nurse = nurses_available[1]
        nurses_workload[nurse, d] = sum(room_workload[:, d])
        for r in 1:R
          assignment_matrix[nurse, r, d] = true
          println("Nurse $nurse has been assigned to room $r at day $d")
        end
      else
        for n in nurses_available
          personal_workload = 0
          for r in 1:R
            if any(assignment_matrix[:, r, d]) #room already taken
              continue
            end
            if personal_workload + room_workload[r, d] <= (nurse_available_on_day_load[n])[d-1] && skill_level_by_room[r,d] <= nurse_skill_level[n]
              personal_workload += room_workload[r, d] 
              assignment_matrix[n, r, d] = true
              println("Nurse $n has been assigned to room $r at day $d")
            end
          end
          nurses_workload[n, d] = personal_workload
        end
      end
    end

    #if there are still rooms without nurses while respecting workload but no more skill level
    not_allocated_rooms_day = [(r, d) for r in 1:R, d in 1:D if !any(assignment_matrix[:, r, d])]
    println(not_allocated_rooms_day)
    if length(not_allocated_rooms_day) > 0
      days_with_remaining_rooms = unique([j[2] for j in not_allocated_rooms_day])
      for d in days_with_remaining_rooms
        nurses_available = [n for n in 1:N if haskey(nurse_available_on_day_load[n], d - 1)]
        rooms_list = [j[1] for j in not_allocated_rooms_day if j[2] == d]
        for n in nurses_available
          personal_workload = nurses_workload[n, d]
          for r in rooms_list
            if personal_workload + room_workload[r, d] <= (nurse_available_on_day_load[n])[d-1]
              personal_workload += room_workload[r, d]
              assignment_matrix[n, r, d] = true
              println("Nurse $n has been assigned to room $r at day $d")
            end
          end
          nurses_workload[n, d] = personal_workload
        end
      end
    end

    #if there are still rooms without nurses no attention on skill level nor workload
    not_allocated_rooms_day = [(r, d) for r in 1:R, d in 1:D if !any(assignment_matrix[:, r, d])]
    println(not_allocated_rooms_day)
    if length(not_allocated_rooms_day) > 0
      days_with_remaining_rooms = unique([j[2] for j in not_allocated_rooms_day])
      for d in days_with_remaining_rooms
        nurses_available = [n for n in 1:N if haskey(nurse_available_on_day_load[n], d - 1)]
        rooms_list = [j[1] for j in not_allocated_rooms_day if j[2] == d]
        for n in nurses_available
          personal_workload = nurses_workload[n, d]
          for r in rooms_list
            personal_workload += room_workload[r, d]
            assignment_matrix[n, r, d] = true
            println("Nurse $n has been assigned to room $r at day $d")
          end
          nurses_workload[n, d] = personal_workload
        end
      end
    end
end 

Heursitic("test02.json")
