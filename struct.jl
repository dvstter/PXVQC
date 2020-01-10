export MPS, QuantumCircuit, OnebodyGate, TwobodyGate, Variable, OnebodyVGate, TwobodyVGate

#=
    Quantum Structures.
    For conventions, please refer to the original implementation
=#
struct MPS
    data::Vector{Array{Number, 3}}
end

Base.length(sys::MPS) = return length(sys.data)

abstract type AbstractGate end
abstract type AbstractOnebodyGate <: AbstractGate end
abstract type AbstractTwobodyGate <: AbstractGate end

struct OnebodyGate <: AbstractOnebodyGate
    index::Int64
    operation::Array{Number, 2}
    description::String
    OnebodyGate(idx, op, gatename) = size(op) != (2,2) ? error("OnebodyGate parameter wrong") : new(idx, op, "$gatename -> $idx")
end

struct TwobodyGate <: AbstractTwobodyGate
    indexes::Array{Int64, 1}
    operation::Array{Number, 4}
    description::String
    TwobodyGate(idxs, op, gatename) = size(op) != (2,2,2,2) || length(idxs) != 2 ? error("TwobodyGate parameter wrong") : new(idxs, op, "$gatename $(idxs[1]) -> $(idxs[2])")
end

mutable struct Variable
    var::Number
end

struct OnebodyVGate <: AbstractOnebodyGate
    index::Int64
    variable::Variable
    op_func::Any
    description::String
    OnebodyVGate(idx, var, op_f, gatename) = new(idx, var, op_f, "(v)$gatename($var) -> idx")
end

struct TwobodyVGate <: AbstractTwobodyGate
    indexes::Array{Int64, 1}
    variable::Variable
    op_func::Any
    description::String
    TwobodyVGate(idxs, var, op_f, gatename) = new(idxs, var, op_f, "(v)$gatename($var) $(idxs[1]) -> $(idxs[2])")
end

struct QuantumCircuit
    gates::Array{AbstractGate, 1}
end

# important overloading method for QuantumCircuit
Base.iterate(circuit::QuantumCircuit, index::Int64) = index > length(circuit.gates) ? nothing : (circuit.gates[index], index+1)
Base.length(circuit::QuantumCircuit, gatetype::Symbol=:all) = begin
    len = length(circuit.gates)
    if gatetype == :variable
        len = length(filter(g -> typeof(g) <: Union{OnebodyVGate, TwobodyVGate}, circuit.gates))
    end
    return len
end
Base.eltype(::Type{QuantumCircuit}) = AbstractGate
Base.getindex(circuit::QuantumCircuit, index::Int64) = begin
    1 <= index <= length(circuit.gates) || throw(BoundsError(circuit, index))
    return circuit.gates[index]
end
Base.firstindex(circuit::QuantumCircuit) = 1
Base.lastindex(circuit::QuantumCircuit) = length(circuit.gates)
Base.setindex!(circuit::QuantumCircuit, gate::AbstractGate, index::Int64) = begin
    1 <= index <= length(circuit.gates) || throw(BoundsError(circuit, index))
    old_gate = circuit.gates[index]
    circuit.gates[index] = gate
    return old_gate
end
Base.push!(circuit::QuantumCircuit, gate::AbstractGate) = begin push!(circuit.gates, gate) end
Base.show(io::IO, circuit::QuantumCircuit) = begin
    for each in circuit
        print(io, each.description)
    end
end
