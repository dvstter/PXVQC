export initmps, initcircuit, apply!, 𝑿!, 𝒀!, 𝒁!, 𝑯!, projection, projection_all, measure!

function initmps(basis::AbstractArray)
    res = []
    for x in basis
        t = []
        if x == 0
            t = reshape([1,0], 1, 2, 1)
        elseif x == 1
            t = reshape([0,1], 1, 2, 1)
        else
            x = rand()
            y = sqrt(1-x)
            t = reshape([sqrt(x), y], 1, 2, 1)
        end
        push!(res, Number.(t))
    end
    return MPS(res)
end

function initcircuit()
    return QuantumCircuit([])
end

# get one controlled single-qubit gate
# note that this gate can only be applied on two adjacent site
function __reverse_gate(gate::AbstractArray)
    size(gate) != (2,2,2,2) && error("parameter wrong")
    return permutedims(gate, [2,1,4,3])
end

# get one swap gate
# note that this gate can only be applied on two adjacent site
function __swap_gate()
    cnot₁ = __controlled_gate([0 1;1 0])
    cnot₂ = deepcopy(cnot₁)
    cnot₃ = deepcopy(cnot₁)
    g₂ = *(cnot₁, cnot₂, (3,4), (2,1))
    return *(g₂, cnot₃, (4,3), (1,2))
end

# get one controlled single-qubit gate
# note that this gate can only be applied on two adjacent site
function __controlled_gate(gate::AbstractArray)
    size(gate) !== (2,2) && error("parameter wrong")
    # generate |0><0|⨂Identity
    a = *([1 0;0 0], [1 0;0 1], (), ())
    # generate |1><1|⨂OnebodyGate
    b = *([0 0;0 1], gate, (), ())
    return permutedims(a+b, [1,3,2,4])
end

function __adjacent_SWAP(sys::MPS, i₁::Int64, i₂::Int64)
    abs(i₂ - i₁) != 1 && error("parameter wrong")
    i₁ > i₂ && begin i₁, i₂ = i₂, i₁ end
    g₂ = __swap_gate()
    sites = *(sys.data[i₁], sys.data[i₂], (3,), (1,))
    sites = *(sites, g₂, (2,3), (3,4))
    u, vt = tsvd(sites, (4,2))
    sys.data[i₁] = u
    sys.data[i₂] = vt
    return (u, vt)
end

function __apply_onebody_gate(sys::MPS, gate::AbstractOnebodyGate)
    operation = []
    if isa(gate, OnebodyVGate) operation = gate.op_func(gate.variable.var) else operation = gate.operation end
    sys.data[gate.index] = permutedims(*(sys.data[gate.index], operation, (2,), (2,)), [1,3,2])
end

function __apply_twobody_gate(sys::MPS, gate::AbstractTwobodyGate)
    idx₁ = gate.indexes[1]
    idx₂ = gate.indexes[2]
    g₂ = []
    if isa(gate, TwobodyVGate) g₂ = gate.op_func(gate.variable.var) else g₂ = gate.operation end
    idx₁ > idx₂ && begin idx₁, idx₂ = idx₂, idx₁; g₂ = __reverse_gate(g₂) end
    sorder = [(x, x+1) for x in idx₁:idx₂-2]
    revorder = reverse(sorder)
    for (x, y) in sorder
        __adjacent_SWAP(sys, x, y)
    end
    sites = *(sys.data[idx₂-1], sys.data[idx₂], (3,), (1,))
    sites = *(sites, g₂, (2,3), (3,4))
    u, vt = tsvd(sites, (4,2))
    sys.data[idx₂-1] = u
    sys.data[idx₂] = vt
    for (x, y) in revorder
        __adjacent_SWAP(sys, x, y)
    end
end

function apply!(sys::MPS, circuit::QuantumCircuit, range::Union{Nothing, UnitRange}=nothing)
    gates = circuit.gates
    if range != nothing gates = gates[range] end
    for x in gates
        isa(x, OnebodyGate) && __apply_onebody_gate(sys, x)
        isa(x, TwobodyGate) && __apply_twobody_gate(sys, x)
    end
    return sys
end

# operators overloading
# |> operator will apply on the mps itself, but * will generate one new mps state(which means not affect the original mps system)
function Base.:|>(circuit::QuantumCircuit, sys::MPS)
    apply!(sys, circuit)
end

function Base.:*(sys::MPS, circuit::QuantumCircuit)
    tempsys = deepcopy(sys)
    return apply!(tempsys, circuit)
end

function Base.:*(circuit::QuantumCircuit, sys::MPS)
    tempsys = deepcopy(sys)
    return apply!(tempsys, circuit)
end

# directly apply gate on the mps
𝑿!(sys::MPS, index::Int64) = __apply_onebody_gate(sys, pauli_x(index))
𝒀!(sys::MPS, index::Int64) = __apply_onebody_gate(sys, pauli_y(index))
𝒁!(sys::MPS, index::Int64) = __apply_onebody_gate(sys, pauli_z(index))
𝑯!(sys::MPS, index::Int64) = __apply_onebody_gate(sys, hadamard(index))

function projection(sys::MPS, basis::AbstractArray, output::Bool=false)
    length(basis) != length(sys.data) && error("parameter wrong, projection and sys's dims not equal")
    proj = initmps(basis)
    nqubits = length(sys.data)
    t = *(sys.data[nqubits], proj.data[nqubits], (2,), (2,))
    for x in (nqubits-1):-1:1
        t₀ = *(sys.data[x], proj.data[x], (2,), (2,))
        t = *(t, t₀, (1,3), (2,4))
        # these three methods should be contract to one method, but i don't know how to yet
        t = permutedims(t, (2,1,3,4))
        t = permutedims(t, (1,2,4,3))
        t = permutedims(t, (4,2,3,1))
    end
    res = reshape(t, 1)[1]
    output && println("projection for |ϕ>=|", join(basis, ""), "> == ", res)
    return res
end

function projection_all(sys::MPS, message::String="")
    println(message)
    nqubits = length(sys.data)
    for x in 0:(2^nqubits - 1)
        n = reverse(digits(x, base=2))
        basis = [[0 for _ in 1:(nqubits-length(n))];n]
        projection(sys, basis, true)
    end
end

function measure!(sys::MPS, i::Int64, output_probability::Bool=false)
    # restore collapsed state for |0> and |1>
    state_collapse₀ = permutedims(*(sys.data[i], [1 0;0 0], (2,), (2,)), [1,3,2])
    state_collapse₁ = permutedims(*(sys.data[i], [0 0;0 1], (2,), (2,)), [1,3,2])

    # assume collapse |0> then get that probability
    sys.data[i] = state_collapse₀
    nqubits = length(sys.data)
    t = *(sys.data[nqubits], sys.data[nqubits], (2,), (2,))
    for x in (nqubits-1):-1:1
        t₀ = *(sys.data[x], sys.data[x], (2,), (2,))
        t = *(t, t₀, (1,3), (2,4))
        # these three methods should be contract to one method, but i don't know how to yet
        t = permutedims(t, (2,1,3,4))
        t = permutedims(t, (1,2,4,3))
        t = permutedims(t, (4,2,3,1))
    end

    probability₀ = reshape(t, 1)[1]
    if output_probability
        println("probability for |0> : ", probability₀)
    end

    # collapse #i qubit and return the classical information
    if rand() < probability₀
        sys.data[i] = state_collapse₀./sqrt(probability₀)
        return 0
    else
        sys.data[i] = state_collapse₁./sqrt(1-probability₀)
        return 1
    end
end
