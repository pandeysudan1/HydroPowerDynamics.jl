# Trollheim-inspired TurbineLookup step test
# Gate/load command: 0.5 -> 0.7 pu at t = 5 s
#
# The turbine lookup itself is algebraic. Head and shaft speed are fixed here
# so this example isolates the static hill-chart response. Governor/servo
# dynamics are intentionally excluded and will be added later.

using HydroPowerDynamics
using ModelingToolkit
using OrdinaryDiffEq
using Printf

const RHO = 1000.0
const G = 9.81
const H = 371.0
const D = 2.5
const N11 = 50.0
const N_RPM = N11 * sqrt(H) / D
const OMEGA = 2pi * N_RPM / 60
const TSTEP = 5.0
const GATE0 = 0.5
const GATE1 = 0.7

@mtkmodel HydraulicBoundaryZ begin
    @parameters begin
        p = 101_325.0
        z = 0.0
    end
    @components port = HydraulicPortZ()
    @equations begin
        port.p ~ p
        port.z ~ z
    end
end

@mtkmodel ShaftSpeedBoundary begin
    @parameters omega0 = OMEGA
    @components port = RotationalPort()
    @equations begin
        port.omega ~ omega0
        port.phi ~ omega0 * t
    end
end

@mtkmodel GateStepSource begin
    @parameters begin
        gate0 = GATE0
        gate1 = GATE1
        tstep = TSTEP
    end
    @variables x(t) = 0.0
    @components out = SignalOutPort()
    @equations begin
        D(x) ~ 0.0
        out.u ~ ifelse(t < tstep, gate0, gate1)
    end
end

@mtkcompile sys begin
    upper = HydraulicBoundaryZ(p = 101_325.0 + RHO*G*H, z = 0.0)
    turbine = TurbineLookup(rho=RHO, g=G, D_runner=D)
    lower = HydraulicBoundaryZ(p = 101_325.0, z = 0.0)
    shaft = ShaftSpeedBoundary()
    gate = GateStepSource()

    @equations begin
        connect(upper.port, turbine.port_a)
        connect(turbine.port_b, lower.port)
        connect(turbine.shaft, shaft.port)
        connect(gate.out, turbine.gate)
    end
end

report = model_structure_report(sys)
println("STRUCTURE=", report)

prob = ODEProblem(sys, [sys.gate.x => 0.0], (0.0, 10.0))
sol = solve(prob, Rodas5P(); tstops=[TSTEP], saveat=0.1,
            reltol=1e-10, abstol=1e-12)

@printf("Pre-step:  Q=%.6f m3/s eta=%.6f Pm=%.6f MW\n",
    sol[sys.turbine.Q][1],
    sol[sys.turbine.eta][1],
    sol[sys.turbine.P_mech][1]/1e6)

@printf("Post-step: Q=%.6f m3/s eta=%.6f Pm=%.6f MW\n",
    sol[sys.turbine.Q][end],
    sol[sys.turbine.eta][end],
    sol[sys.turbine.P_mech][end]/1e6)
