export collect_parameters, putback_parameters!

function collect_parameters(circuit::QuantumCircuit)::Array{Number, 1}
    res  = []
    for gate in circuit.gates
        if isa(gate, OnebodyVGate) || isa(gate, TwobodyVGate)
            push!(res, gate.variable.var)
        end
    end
    return res
end

function putback_parameters!(circuit::QuantumCircuit, vars::AbstractArray)
    idx = 1
    for gate in circuit.gates
        if isa(gate, OnebodyVGate) || isa(gate, TwobodyVGate)
            gate.variable.var = vars[idx]
            idx+=1
        end
    end
end
#=
@adjoint *(sys::MPS, circuit::QuantumCircuit) = begin
    r = *(sys, circuit)

end=#
