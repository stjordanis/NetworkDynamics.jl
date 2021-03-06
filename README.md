# NetworkDynamics

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://fhell.github.io/NetworkDynamics.jl/dev)
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://fhell.github.io/NetworkDynamics.jl/stable) -->

***
News: We just released a preprint paper about NetworkDynamics.jl on [arXiv](https://arxiv.org/abs/2012.12696). Make sure to use ND v0.5 if you want to reproduce the results shown there. For more details see also the latest docs ([dev](https://fhell.github.io/NetworkDynamics.jl/dev)).
***

A package for working with dynamical systems on complex networks. NetworkDynamics.jl provides an interface between [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). It allows to define several types of dynamic and static nodes and edges and to link them up in order to create complex network dynamics.

The behavior of a node or an edge can be described by algebraic equations, by differential algebraic equation (DAEs) in mass matrix form, by ordinary differential equations (ODE) or by delay differential equations (DDE). Stochastic ordinary differential equations (SDE) can be implemented as a [two-layer network](https://github.com/FHell/NetworkDynamics.jl/blob/master/examples/sde.jl). For details see the [docs](https://fhell.github.io/NetworkDynamics.jl/dev).

## Installation

Installation is straightforward with Julia's package manager.

```julia-repl
(v1.5) pkg> add NetworkDynamics
```

## Getting started

Check out our step-by-step tutorial as a [jupyter notebook](https://github.com/FHell/NetworkDynamics.jl/blob/master/examples/getting_started_with_network_dynamics.ipynb) or [in the docs](https://fhell.github.io/NetworkDynamics.jl/dev/getting_started_with_network_dynamics/).

## PowerDynamics

[PowerDynamics.jl](https://juliaenergy.github.io/PowerDynamics.jl/stable/) is an open-source framework for dynamic power grid modeling and analysis build on top of NetworkDynamics.jl.
