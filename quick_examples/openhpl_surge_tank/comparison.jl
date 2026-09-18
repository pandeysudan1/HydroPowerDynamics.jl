# Reduced vs OpenHPL-style surge-tank comparison
#
# This example compares the constitutive structure of the existing one-state
# SurgeTank with the new two-state OpenHPLSurgeTank.
#
# The reduced model maps flow directly to level and pressure.
# The OpenHPL-style model retains water-column momentum, so bottom pressure
# and shaft flow interact dynamically.

using HydroPowerDynamics
using ModelingToolkit
using OrdinaryDiffEq
using Plots

# A pressure boundary that changes at t = 5 s.
@mtkmodel PressureBoundary begin
    @parameters begin
        p0 = 101_325.0 + 1000.0 * 9.81 * 50.0
        dp = 15_000.0
        t_step = 5.0
    end
    @components begin
        port = HydraulicPort()
    end
    @equations begin
        port.p ~ p0 + dp * ifelse(t >= t_step, 1.0, 0.0)
    end
end

@mtkcompile sys begin
    boundary = PressureBoundary()
    tank = OpenHPLSurgeTank(
        H = 87.0,
        L = 87.0,
        D_tank = 3.4,
        h_0 = 50.0,
        Vdot_0 = 0.0,
        f = 0.02,
    )
    @equations begin
        connect(boundary.port, tank.port)
    end
end

u0 = [
    sys.tank.h => 50.0,
    sys.tank.Vdot => 0.0,
]

prob = ODEProblem(sys, u0, (0.0, 60.0))
sol = solve(prob, Rodas5P(); reltol = 1e-8, abstol = 1e-9)

p1 = plot(sol; idxs = [sys.tank.h], xlabel = "Time [s]", ylabel = "h [m]",
          title = "OpenHPL surge-tank water level", legend = false)
p2 = plot(sol; idxs = [sys.tank.Vdot], xlabel = "Time [s]", ylabel = "Vdot [m^3/s]",
          title = "Water-column momentum response", legend = false)

display(p1)
display(p2)
