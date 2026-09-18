using Test
using HydroPowerDynamics

@testset "Trollheim TurbineLookup 0.5 -> 0.7 gate step" begin
    H = 371.0
    D = 2.5
    rho = 1000.0
    g = 9.81
    n11 = 50.0

    q11_05 = bilinear3x3_clamped(
        0.5, n11,
        0.2,0.6,1.0, 40.0,50.0,60.0,
        0.060,0.062,0.060,
        0.185,0.190,0.184,
        0.300,0.30735108593749,0.298
    )
    eta_05 = bilinear3x3_clamped(
        0.5, n11,
        0.2,0.6,1.0, 40.0,50.0,60.0,
        0.72,0.76,0.70,
        0.89,0.93,0.88,
        0.94,0.97,0.92
    )

    q11_07 = bilinear3x3_clamped(
        0.7, n11,
        0.2,0.6,1.0, 40.0,50.0,60.0,
        0.060,0.062,0.060,
        0.185,0.190,0.184,
        0.300,0.30735108593749,0.298
    )
    eta_07 = bilinear3x3_clamped(
        0.7, n11,
        0.2,0.6,1.0, 40.0,50.0,60.0,
        0.72,0.76,0.70,
        0.89,0.93,0.88,
        0.94,0.97,0.92
    )

    Q05 = q11_05 * D^2 * sqrt(H)
    Q07 = q11_07 * D^2 * sqrt(H)
    Pm05 = eta_05 * rho * g * Q05 * H
    Pm07 = eta_07 * rho * g * Q07 * H

    @test isapprox(q11_05, 0.158; atol=1e-12)
    @test isapprox(eta_05, 0.8875; atol=1e-12)
    @test isapprox(Q05, 19.020593280705; rtol=1e-12)
    @test isapprox(Pm05/1e6, 61.437755012815; rtol=1e-12)

    @test isapprox(q11_07, 0.2193377714843725; atol=1e-12)
    @test isapprox(eta_07, 0.94; atol=1e-12)
    @test isapprox(Q07, 26.4046490031675; rtol=1e-12)
    @test isapprox(Pm07/1e6, 90.333985047907; rtol=1e-12)

    @test Q07 > Q05
    @test Pm07 > Pm05
end
