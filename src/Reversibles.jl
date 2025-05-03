module Reversibles
import FunctionWrappers: FunctionWrapper

export Trail, savestate!, restorestate!
export ReversibleAssignable, update!, registerOnUpdate!, inc!, dec!
export ReversibleInt, ReversibleFloat64, ReversibleString
export ReversibleSet, push!, delete!, registerOnPush!, registerOnDelete!, in, length, iterate
export ReversibleCountingSet, getamount
export show

using DataStructures
import DataStructures: update!, push!, delete!, top

mutable struct Trail 
    const s::Stack{Vector{Function}}
    magic::Int

    function Trail()
        obj = new(Stack{Vector{Function}}(), 0)
        push!(obj.s, Vector{Function}())
        obj
    end
end

function savestate!(trail::Trail)
    trail.magic += 1
    push!(trail.s, Vector{Function}())
end

function restorestate!(trail::Trail)
    tops = pop!(trail.s)
    trail.magic += 1
    for f in Iterators.reverse(tops)
        f()
    end
end

mutable struct ReversibleAssignable{T}
    value::T
    trail::Trail
    lastmagic::Int
    onupdate::Vector{FunctionWrapper{Nothing, Tuple{T, T}}}

    function ReversibleAssignable{T}(value::T, trail::Trail) where {T}
        obj = new{T}(value, trail)
        obj.lastmagic = trail.magic-1
        obj.onupdate = Vector{Function}()
        obj
    end
end

function Base.show(io::IO,x::ReversibleAssignable{T}) where {T}
    print(io, x.value)
end

function registerOnUpdate!(r::ReversibleAssignable{T}, fun::Function)::Nothing where T
    push!(r.onupdate, FunctionWrapper{Nothing, Tuple{T, T}}(fun))
    return
end

function update!(r::ReversibleAssignable{T}, newvalue::T) where {T}
    oldvalue = r.value
    if oldvalue == newvalue
        return
    end

    if r.lastmagic != r.trail.magic
        push!(top(r.trail.s), () -> r.value = oldvalue)
        r.lastmagic = r.trail.magic
    end

    r.value = newvalue
    for f in r.onupdate
        f(oldvalue, newvalue)
    end
end

function inc!(r::ReversibleAssignable{T}, x::T) where {T}
    update!(r, r.value + x)
    return r
end

function dec!(r::ReversibleAssignable{T}, x::T) where {T}
    update!(r, r.value - x)
    return r
end

ReversibleInt = ReversibleAssignable{Int}
ReversibleFloat64 = ReversibleAssignable{Float64}
ReversibleString = ReversibleAssignable{String}

mutable struct ReversibleSet{T}
    value::Set{T}
    trail::Trail
    lastmagic::Int
    onpush::Vector{FunctionWrapper{Nothing, Tuple{Set{T}, T}}}
    ondelete::Vector{FunctionWrapper{Nothing, Tuple{Set{T}, T}}}

    function ReversibleSet{T}(set::Set{T}, trail::Trail) where {T}
        new(set, trail, trail.magic, Vector{FunctionWrapper{Nothing, Tuple{Set{T}, T}}}(), Vector{FunctionWrapper{Nothing, Tuple{Set{T}, T}}}())
    end

    function ReversibleSet{T}(trail::Trail) where {T}
        ReversibleSet{T}(Set{T}(), trail)
    end
end

function registerOnPush!(r::ReversibleSet{T}, fun::Function) where T
    push!(r.onpush, FunctionWrapper{Nothing, Tuple{Set{T}, T}}(fun))
end

function registerOnDelete!(r::ReversibleSet{T}, fun::Function) where T
    push!(r.ondelete, FunctionWrapper{Nothing, Tuple{Set{T}, T}}(fun))
end

function Base.in(r::ReversibleSet{T}, value::T) where {T}
    value in r.value
end

function Base.length(r::ReversibleSet{T}) where {T}
    length(r.value)
end

function Base.iterate(r::ReversibleSet{T}) where {T}
    iterate(r.value)
end

function Base.iterate(r::ReversibleSet{T}, state) where {T}
    iterate(r.value, state)
end

function push!(r::ReversibleSet{T}, value::T) where {T}
    if value in r.value
        return
    end

    if r.lastmagic != r.trail.magic
        cvalue = copy(r.value)
        push!(top(r.trail.s), () -> r.value = cvalue)
        r.lastmagic = r.trail.magic
    end
    push!(r.value, value)

    for f in r.onpush
        f(r.value, value)
    end
end

function delete!(r::ReversibleSet{T}, value::T) where {T}
    if value ∉ r.value
        return
    end

    if r.lastmagic != r.trail.magic
        cvalue = copy(r.value)
        push!(top(r.trail.s), () -> r.value = cvalue)
        r.lastmagic = r.trail.magic
    end
    pop!(r.value, value)

    for f in r.ondelete
        f(r.value, value)
    end

    return r
end

function Base.show(io::IO,x::ReversibleSet{T}) where {T}
    print(io, x.value)
end

mutable struct ReversibleCountingSet{T}
    value::Dict{T,Int}
    trail::Trail
    lastmagic::Int
    onpush::Vector{FunctionWrapper{Nothing, Tuple{Dict{T,Int}, T}}}
    ondelete::Vector{FunctionWrapper{Nothing, Tuple{Dict{T,Int}, T}}}

    function ReversibleCountingSet{T}(init::Dict{T,Int}, trail::Trail) where {T}
        new(init, trail, trail.magic, Vector{FunctionWrapper{Nothing, Tuple{Dict{T,Int}, T}}}(), Vector{FunctionWrapper{Nothing, Tuple{Dict{T,Int}, T}}}())
    end

    function ReversibleCountingSet{T}(trail::Trail) where {T}
        ReversibleCountingSet{T}(Dict{T,Int}(), trail)
    end
end

function registerOnPush!(r::ReversibleCountingSet{T}, fun::Function) where T
    push!(r.onpush, FunctionWrapper{Nothing, Tuple{Dict{T,Int}, T}}(fun))
end

function registerOnDelete!(r::ReversibleCountingSet{T}, fun::Function) where T
    push!(r.ondelete, FunctionWrapper{Nothing, Tuple{Dict{T,Int}, T}}(fun))
end

function Base.in(r::ReversibleCountingSet{T}, value::T) where {T}
    value in r.value
end

function Base.length(r::ReversibleCountingSet{T}) where {T}
    length(r.value)
end

function Base.iterate(r::ReversibleCountingSet{T}) where {T}
    iterate(keys(r.value))
end

function Base.iterate(r::ReversibleCountingSet{T}, state) where {T}
    iterate(keys(r.value), state)
end

function increment_count!(dict::Dict{T, Int}, key::T) where T
    dict[key] = get!(dict, key, 0) + 1
end

function decrement_count!(dict::Dict{T, Int}, key::T) where T
    if haskey(dict, key)
        dict[key] -= 1
        if dict[key] == 0
            delete!(dict, key)
        end
    end
end

function push!(r::ReversibleCountingSet{T}, value::T) where {T}
    if r.lastmagic != r.trail.magic
        cvalue = copy(r.value)
        push!(top(r.trail.s), () -> r.value = cvalue)
        r.lastmagic = r.trail.magic
    end
    increment_count!(r.value, value)

    if r.value[value] == 1
        for f in r.onpush
            f(r.value, value)
        end
    end
end

function delete!(r::ReversibleCountingSet{T}, value::T) where {T}
    if value ∉ keys(r.value)
        return
    end

    if r.lastmagic != r.trail.magic
        cvalue = copy(r.value)
        push!(top(r.trail.s), () -> r.value = cvalue)
        r.lastmagic = r.trail.magic
    end
    decrement_count!(r.value, value)

    if value ∉ keys(r.value)
        for f in r.ondelete
            f(r.value, value)
        end
    end

    return r
end

function getamount(r::ReversibleCountingSet{T}, value::T) where T
    if value ∉ keys(r.value)
        return 0
    end
    return r.value[value]
end

function Base.show(io::IO,x::ReversibleCountingSet{T}) where {T}
    print(io, x.value)
end

end