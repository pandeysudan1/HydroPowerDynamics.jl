using Test
using ModelingToolkit
using HydroPowerDynamics

@testset "Physical connector interfaces" begin

    @mtkmodel HydraulicPressureBoundary begin
        @components port = HydraulicPort()
        @equations port.p ~ 2.0e5
    end

    @mtkmodel HydraulicFlowLoad begin
        @components port = HydraulicPort()
        @equations port.dm ~ 10.0
    end

    @mtkcompile hydraulic_sys begin
        source = HydraulicPressureBoundary()
        load = HydraulicFlowLoad()
        @equations connect(source.port, load.port)
    end
    @test hydraulic_sys !== nothing

    @mtkmodel RotationalSpeedBoundary begin
        @components port = RotationalPort()
        @equations begin
            port.phi ~ 0.0
            port.omega ~ 100.0
        end
    end

    @mtkmodel TorqueLoad begin
        @components port = RotationalPort()
        @equations port.tau ~ 25.0
    end

    @mtkcompile rotational_sys begin
        source = RotationalSpeedBoundary()
        load = TorqueLoad()
        @equations connect(source.port, load.port)
    end
    @test rotational_sys !== nothing

    @mtkmodel TranslationalMotionBoundary begin
        @components port = TranslationalPort()
        @equations begin
            port.x ~ 0.0
            port.v ~ 0.0
        end
    end

    @mtkmodel ForceLoad begin
        @components port = TranslationalPort()
        @equations port.f ~ 100.0
    end

    @mtkcompile translational_sys begin
        source = TranslationalMotionBoundary()
        load = ForceLoad()
        @equations connect(source.port, load.port)
    end
    @test translational_sys !== nothing

    @mtkmodel ACVoltageBoundary begin
        @components port = ElectricalACPort()
        @equations begin
            port.vr ~ 400.0
            port.vi ~ 0.0
        end
    end

    @mtkmodel ACCurrentLoad begin
        @components port = ElectricalACPort()
        @equations begin
            port.ir ~ 10.0
            port.ii ~ -2.0
        end
    end

    @mtkcompile ac_sys begin
        source = ACVoltageBoundary()
        load = ACCurrentLoad()
        @equations connect(source.port, load.port)
    end
    @test ac_sys !== nothing

    @mtkmodel DCVoltageBoundary begin
        @components port = ElectricalDCPort()
        @equations port.v ~ 800.0
    end

    @mtkmodel DCCurrentLoad begin
        @components port = ElectricalDCPort()
        @equations port.i ~ 5.0
    end

    @mtkcompile dc_sys begin
        source = DCVoltageBoundary()
        load = DCCurrentLoad()
        @equations connect(source.port, load.port)
    end
    @test dc_sys !== nothing
end
