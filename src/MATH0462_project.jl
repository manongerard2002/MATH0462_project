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
export LSManager, getNewSolution!, addSolution!, addWaitingFeasibleSolution!, getWaitingFeasibleSolution!, getBestSolution

include("local_search/LNSNurses.jl")
include("local_search/PatientOPT.jl")
export LNSNurses, PatientOPT

include("validator.jl")
export get_violations_scores

end