export pauli_x, pauli_y, pauli_z, hadamard, cnot, ry, ∂ry

# frequently-used quantum gates
pauli_x(i::Int64) = OnebodyGate(i, [0 1;1 0])
pauli_y(i::Int64) = OnebodyGate(i, [0 -im;im 0])
pauli_z(i::Int64) = OnebodyGate(i, [1 0;0 -1])
hadamard(i::Int64) = OnebodyGate(i, 1/sqrt(2)*[1 1;1 -1])
cnot(con::Int64, op::Int64) = TwobodyGate([con, op], __controlled_gate([0 1;1 0]))
ry(i::Int64, θ::Union{Number, Variable}) = begin
    if isa(θ, Variable) OnebodyVGate(i, θ, 𝒗 -> [cos(𝒗) -sin(𝒗);sin(𝒗) cos(𝒗)])
    else OnebodyGate(i, [cos(θ) -sin(θ);sin(θ) cos(θ)]) end
end
∂ry(i::Int64, θ::Union{Number, Variable}) = begin
    if isa(θ, Variable) return OnebodyVGate(i, θ, 𝒗 -> [cos(𝒗+1/2pi) -sin(𝒗+1/2pi);sin(𝒗+1/2pi) cos(𝒗+1/2pi)]) end
    return ry(i, θ)
end
