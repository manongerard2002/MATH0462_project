module PatientOPT
using ..Reversibles
using ..LocalSearchUtils
using Combinatorics

function patient_2_opt(startTime, limitTime, m, ordering::LSOrdering; until_convergence=true)
    patients_g1 = [p for p in ordering.patients if m.patient_gender[p] == 1]
    patients_g2 = [p for p in ordering.patients if m.patient_gender[p] == 2]
    patients_alone = () -> [p for p in ordering.patients if m.v_patient_admission_day[p].value != m.D + 1 && m.c_patients_present_in_room_n[m.v_patient_room[p].value].value == 1]
    if until_convergence
        return run_until_convergence(startTime, limitTime, [
            () -> patient_2_opt_p(startTime, limitTime, m, combinations(patients_g1, 2)),
            () -> patient_2_opt_p(startTime, limitTime, m, combinations(patients_g2, 2)),
            () -> patient_2_opt_p(startTime, limitTime, m, combinations(patients_alone(), 2))
        ])
    else
        m1 = patient_2_opt_p(startTime, limitTime, m, combinations(patients_g1, 2), until_convergence=until_convergence)
        m2 = patient_2_opt_p(startTime, limitTime, m, combinations(patients_g2, 2), until_convergence=until_convergence)
        m3 = patient_2_opt_p(startTime, limitTime, m, combinations(patients_alone(), 2), until_convergence=until_convergence)
        return m1 || m2 || m3
    end
end

function exchange_patients(m, p1, p2)
    if p1 == p2
        return false
    end
    (d1, d2) = m.v_patient_admission_day[p1].value, m.v_patient_admission_day[p2].value
    (r1, r2) = m.v_patient_room[p1].value, m.v_patient_room[p2].value

    orig_obj = m.s_total.value
    orig_vio = m.e_total.value

    update!(m.v_patient_admission_day[p1], d2)
    update!(m.v_patient_admission_day[p2], d1)
    update!(m.v_patient_room[p1], r2)
    update!(m.v_patient_room[p2], r1)

    if is_improving(m, orig_vio, orig_obj)
        return true
    else
        # revert state
        update!(m.v_patient_admission_day[p1], d1)
        update!(m.v_patient_admission_day[p2], d2)
        update!(m.v_patient_room[p1], r1)
        update!(m.v_patient_room[p2], r2)
        return false
    end
end

function patient_2_opt_p(startTime, limitTime, m, pairs; until_convergence=true)
    if until_convergence
        return run_until_convergence(startTime, limitTime, [
            () -> exchange_patients(m, p1, p2) for (p1, p2) in pairs
        ])
    else
        modified = false
        for (p1, p2) in pairs
            modified |= exchange_patients(m, p1, p2)
        end
        return modified
    end
end

function exchange_patients_n_opt(m, patients...)
    ds = [m.v_patient_admission_day[p].value for p in patients]
    rs = [m.v_patient_room[p].value for p in patients]

    orig_obj = m.s_total.value
    orig_vio = m.e_total.value

    for i in eachindex(patients)
        next_i = (i % length(patients)) + 1
        update!(m.v_patient_admission_day[patients[i]], ds[next_i])
        update!(m.v_patient_room[patients[i]], rs[next_i])
    end

    if is_improving(m, orig_vio, orig_obj)
        return true
    else
        # revert state
        for i in eachindex(patients)
            update!(m.v_patient_admission_day[patients[i]], ds[i])
            update!(m.v_patient_room[patients[i]], rs[i])
        end
        return false
    end
end

function patient_3_opt(startTime, limitTime, m, ordering::LSOrdering; until_convergence=true)
    patients_g1 = [p for p in ordering.patients if m.patient_gender[p] == 1]
    patients_g2 = [p for p in ordering.patients if m.patient_gender[p] == 2]

    if until_convergence
        return run_until_convergence(startTime, limitTime, [
            () -> patient_3_opt_p(startTime, limitTime, m, combinations(patients_g1, 3)),
            () -> patient_3_opt_p(startTime, limitTime, m, combinations(patients_g2, 3))
        ])
    else
        m1 = patient_3_opt_p(startTime, limitTime, m, combinations(patients_g1, 3), until_convergence=until_convergence)
        m2 = patient_3_opt_p(startTime, limitTime, m, combinations(patients_g2, 3), until_convergence=until_convergence)
        return m1 || m2
    end
end

function patient_3_opt_p(startTime, limitTime, m, triplets; until_convergence=true)
    if until_convergence
        return run_until_convergence(startTime, limitTime, [
            () -> exchange_patients_n_opt(m, p1, p2, p3) for (p1, p2, p3) in triplets
        ])
    else
        modified = false
        for (p1, p2, p3) in triplets
            modified |= exchange_patients_n_opt(m, p1, p2, p3)
        end
        return modified
    end
end

export patient_2_opt, patient_2_opt_p, exchange_patients, exchange_patients_n_opt, patient_3_opt, patient_3_opt_p
end