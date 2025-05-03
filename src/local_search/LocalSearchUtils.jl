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

struct LSManagerListener
    newSolutions::Vector{Model}
    lock::ReentrantLock

    function LSManagerListener()
        new(Vector{Model}(), ReentrantLock())
    end
end
export LSManagerListener

function getNewSolution!(listener::LSManagerListener)
    lock(listener.lock)
    try
        if length(listener.newSolutions) == 0
            return nothing
        end
        return pop!(listener.newSolutions)
    finally
        unlock(listener.lock)
    end
end
export getNewSolution!

struct LSManager
    bestSolutions::Vector{Model}
    waitingFeasibleSolutions::Vector{Model}
    lock::ReentrantLock
    nMaxSolutions::Int
    listeners::Vector{LSManagerListener}
    startTime::Float64
    limitTime::Float64

    function LSManager(nMaxSolutions::Int, limitTime)
        new(Vector{Model}(), Vector{Model}(), ReentrantLock(), nMaxSolutions, Vector{LSManagerListener}(), time(), limitTime)
    end
end
export LSManager

function createListener!(manager::LSManager)
    listener = LSManagerListener()
    lock(manager.lock)
    try
        push!(manager.listeners, listener)
    finally
        unlock(manager.lock)
    end
    return listener
end
export createListener!

function addWaitingFeasibleSolution!(manager::LSManager, m::Model, sorted::Bool=true)
    lock(manager.lock)
    try
        #= println("----------------------------------------------------------\n")
        println(@sprintf("New MIP solution: %i (%.2f s)", m.s_total.value, time() - manager.startTime))
        printStats(m)
        flush(stdout) =#
        if sorted #by default feasible solution are find in a decreasing order by Gurobi
            push!(manager.waitingFeasibleSolutions, m)
        else
            idx = findfirst(x -> m.s_total.value > x.s_total.value, manager.waitingFeasibleSolutions)
            if isnothing(idx)
                push!(manager.waitingFeasibleSolutions, m)
            else
                insert!(manager.waitingFeasibleSolutions, idx, m)
            end
        end
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
        if length(manager.bestSolutions) != 0 && solution.s_total.value >= manager.bestSolutions[end].s_total.value
            return
        end

        #= if length(manager.bestSolutions) == 0 || solution.s_total.value < manager.bestSolutions[1].s_total.value
            println("----------------------------------------------------------")
            println(@sprintf("New best solution: %i (%.2f s)", solution.s_total.value, time() - manager.startTime))
            printStats(solution)
            flush(stdout)
        end =#

        solution = Model(solution) #ensure copy
        
        if length(manager.bestSolutions) == manager.nMaxSolutions
            pop!(manager.bestSolutions)
        end
        push!(manager.bestSolutions, solution)
        sort!(manager.bestSolutions, by=s -> s.s_total.value)
        for listener in manager.listeners
            lock(listener.lock)
            try
                push!(listener.newSolutions, solution)
            finally
                unlock(listener.lock)
            end
        end
    finally
        unlock(manager.lock)
    end
end
export addSolution!

function getBestSolution(manager::LSManager)
    lock(manager.lock)
    try
        return manager.bestSolutions[1]
    finally
        unlock(manager.lock)
    end
end
export getBestSolution

end