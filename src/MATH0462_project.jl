module MATH0462_project

include("Reversibles.jl")
export update!, savestate!, restorestate!, getamount

include("MATH0462Model.jl")
using .MATH0462Model
export Model, printStats, loadSolution!, writeSolution!

include("Heuristic_first_available.jl")
export Heuristic_first_available

include("local_search/LocalSearchUtils.jl")
using .LocalSearchUtils
export run_until_convergence, LSObjective, LSOrdering, randomOrdering
export LSManager, LSManagerListener, getNewSolution!, addSolution!, createListener!, addWaitingFeasibleSolution!, getWaitingFeasibleSolution!, getBestSolution


include("UnschedulePatient.jl")
export UnschedulePatient

include("local_search/LNSNurses.jl")
include("local_search/PatientOPT.jl")
export LNSNurses, PatientOPT

include("validator.jl")
export get_violations_scores

end


#= 
include("local_search/LNSNurseRoom.jl")
include("local_search/LNSNurseDay.jl")
include("local_search/LNSPatientRoom.jl")
include("local_search/StupidGreedy.jl")
include("local_search/GreedyPatient.jl")
include("local_search/RandomInit.jl")
include("local_search/UnschedulePatient.jl")
export LNSNurseRoom, LNSNurseDay, LNSNurses, StupidGreedy, LNSPatientRoom, GreedyPatient, Patient2OPT, RandomInit, UnschedulePatient

end =#