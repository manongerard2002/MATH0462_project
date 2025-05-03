module UnschedulePatient
using ..Reversibles
using ..LocalSearchUtils

function unschedule(m, ordering::LSOrdering)
    #println("unschedule optional patients m.e_total=$(m.e_total) m.s_total=$(m.s_total)")
    anything = false
    for p in [p for p in ordering.patients if p in m.optional_patients]
        if m.v_patient_admission_day[p].value == m.D + 1
            continue
        end
        update!(m.v_patient_admission_day[p], m.D + 1)
        anything = true
    end
    return anything
end

end