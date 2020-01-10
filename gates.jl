export pauli_x, pauli_y, pauli_z, hadamard, cnot, ry, ∂ry, crk

# frequently-used quantum gates
pauli_x(i::Int64) = OnebodyGate(i, [0 1;1 0], "pauli-x")
pauli_y(i::Int64) = OnebodyGate(i, [0 -im;im 0], "pauli-y")
pauli_z(i::Int64) = OnebodyGate(i, [1 0;0 -1], "pauli-z")
hadamard(i::Int64) = OnebodyGate(i, 1/sqrt(2)*[1 1;1 -1], "hadamard")
cnot(con::Int64, op::Int64) = TwobodyGate([con, op], __controlled_gate([0 1;1 0]), "cnot")
crk(con::Int64, op::Int64, k::Int64) = TwobodyGate([con, op], __controlled_gate([1 0;0 exp(2pi*im/(2^k))]), "crk")
ry(i::Int64, θ::Union{Number, Variable}) = begin
    if isa(θ, Variable) OnebodyVGate(i, θ, 𝒗 -> [cos(𝒗) -sin(𝒗);sin(𝒗) cos(𝒗)], "ry")
    else OnebodyGate(i, [cos(θ) -sin(θ);sin(θ) cos(θ)], "ry") end
end
