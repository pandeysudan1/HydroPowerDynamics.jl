FROM julia:1.11-bookworm

ENV JULIA_PROJECT=/workspace/HydroPowerDynamics.jl
ENV GKSwstype=100

WORKDIR /workspace/HydroPowerDynamics.jl
COPY . /workspace/HydroPowerDynamics.jl

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

CMD ["bash", "railway/run_agc_trollheim.sh"]
