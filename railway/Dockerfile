FROM julia:1.11-bookworm

ENV JULIA_PROJECT=/workspace/HydroPowerDynamics.jl
ENV JULIA_NUM_PRECOMPILE_TASKS=1
ENV JULIA_PKG_PRECOMPILE_AUTO=0
ENV JULIA_CPU_TARGET=generic
ENV GKSwstype=100

WORKDIR /workspace/HydroPowerDynamics.jl
COPY . /workspace/HydroPowerDynamics.jl

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' \
    && find /root/.julia/compiled -name '*.pidfile' -delete

CMD ["bash", "railway/run_agc_trollheim.sh"]
