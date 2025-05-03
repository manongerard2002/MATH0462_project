module LocalSearchUtils
using Gurobi
using ..MATH0462Model
using Random, Printf

Random.seed!(42)

function is_improving(m, old_vio, old_obj)
    new_vio = m.e_total.value
    new_obj = m.s_total.value
    return new_vio < old_vio || (old_vio == 0 && new_vio == 0 && new_obj < old_obj)
end
export is_improving

function run_until_convergence(startTime, timeLimit, functions; idempotent=true, second_stop_condition=() -> false, after_each=() -> nothing)
    found_something_at_all = false
    i = 1
    n = length(functions)
    time_since_last_update = -1
    limit = idempotent ? n - 1 : n
    while time_since_last_update < limit && !second_stop_condition() && time() - startTime < timeLimit
        found_something = functions[i]()
        after_each()
        found_something_at_all |= found_something
        if found_something
            time_since_last_update = 0
        else
            time_since_last_update += 1
        end

        i = (i % n) + 1
    end
    return found_something_at_all
end
export run_until_convergence

struct LSOrdering
    patients::Vector{Int}
    rooms::Vector{Int}
    nurses::Vector{Int}
    days::Vector{Int}
end

function randomOrdering(m::Model)
    return LSOrdering(
        shuffle(1:m.P),
        shuffle(1:m.R),
        shuffle(1:m.N),
        shuffle(1:m.D)
    )
end

export LSOrdering, randomOrdering

mutable struct LSManager
    bestSolution::Union{Nothing, Model}
    waitingFeasibleSolutions::Vector{Model}
    lock::ReentrantLock
    startTime::Float64
    limitTime::Float64

    function LSManager(limitTime)
        new(nothing, Vector{Model}(), ReentrantLock(), time(), limitTime)
    end
end
export LSManager

function addWaitingFeasibleSolution!(manager::LSManager, m::Model)
    lock(manager.lock)
    try
        push!(manager.waitingFeasibleSolutions, m)
    finally
        unlock(manager.lock)
    end
end
export addWaitingFeasibleSolution!

function getWaitingFeasibleSolution!(manager::LSManager)
    lock(manager.lock)
    try
        if length(manager.waitingFeasibleSolutions) == 0
            return nothing
        end
        return pop!(manager.waitingFeasibleSolutions)
    finally
        unlock(manager.lock)
    end
end
export getWaitingFeasibleSolution!

function addSolution!(manager::LSManager, solution::Model)
    lock(manager.lock)
    try
        if !isnothing(manager.bestSolution) && solution.s_total.value >= manager.bestSolution.s_total.value
            return
        end
        manager.bestSolution = Model(solution) #ensure copy
    finally
        unlock(manager.lock)
    end
end
export addSolution!

function getBestSolution(manager::LSManager)
    lock(manager.lock)
    try
        return manager.bestSolution
    finally
        unlock(manager.lock)
    end
end
export getBestSolution

end