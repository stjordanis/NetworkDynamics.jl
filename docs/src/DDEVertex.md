# Experimental DDE Tutorial

 An `IJulia` [notebook](https://github.com/FHell/NetworkDynamics.jl/tree/master/examples) corresponding to this tutorial will be available on GitHub soon.

#### Topics covered in this tutorial include:
 * constructing a watts_strogatz graph
 * diffusion delay differential equation (DDE)
 * parameter handling with DDE's
 * stiff equations
 * Kuramoto Delay Model

## Network Diffusion

This example explains the use of delay differential equations (DDE's) in NetworkDynamics.jl by modeling a simple diffusion on an undirected ring network with delay.

Let $g$ be a graph with $N$ nodes and adjacency matrix $A$. Let $v = (v_1, \dots, v_n)$ be a vector of (abstract) temperatures or concentrations at each node $i = 1, \dots, N$. The rate of change of state $v_i$ in this artificial example is described by the delay value of $v_i$ with delay $\Tau$ and its difference with its neighbors with coupling strength $\sigma$. We obtain the following ordinary differential equation

```math
\begin{aligned}
\dot v_i(t) = - v_i(t-\Tau) - \sigma * \sum_{i=1}^N A_{ij} (v_i(t) - v_j(t))
\end{aligned}
```

The sum on the right hand side plays the role of a (discrete) gradient. If the temperature at node $i$ is higher than at its neighboring node $j$ it will decrease along that edge.

## Modeling diffusion with delay in NetworkDynamics.jl

From the equation above we see that in this model the dynamics of the nodes consist of internal dynamics as well as a coupling term with the neighboring nodes. In NetworkDynamics.jl the interactions with the neighbors are described by equations for the edges.

```@example DDEVertex
function diffusionedge!(e, v_s, v_d, p, t)
   # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
   e .= .1 * (v_s .- v_d)
   nothing
end
nothing # hide
```
The internal dynamics are determined by the delay value $\dot v_i = v_i(\Tau)$ and are described in the vertex function with help of the history array $h_v$ containing the delay values of the vertex.
```@example DDEVertex
function diffusionvertex!(dv, v, e_s, e_d, h_v, p, t)
   # usually dv, v, e_s, e_d, h_v are arrays, hence we use the broadcasting operator .
   # h_v is the history array of the vertices
   dv .= -h_v
   # edges for which v is the source
   for e in e_s
       dv .-= e
   end
   # edges for which v is the destination
   for e in e_d
       dv .+= e
   end
   nothing
end
nothing # hide
```
# Constructing the Network Dynamics

Constructing the network dynamics is straightforward
```@example DDEVertex
using LightGraphs

### Defining a graph

N = 20 # number of nodes
k = 8  # average degree
g = watts_strogatz(N, k, 0.) # ring-network
```
The Watts-Strogatz algorithm constructs a regular-ring network with $N$ nodes, $k$ neighbors and a probability $p$ of rewiring links, which is chosen as $p=0$ here.

```@example DDEVertex
using NetworkDynamics

nd_diffusion_vertex = DDEVertex(f! = diffusionvertex!, dim = 1)
nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)

nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)
```

Now we hand over the functions we have defined above to the constructors `DDEVertex` and `StaticEdge`, adding information on the dimension of variables at each edge or node with keyword **`dim`**. The wrapper function `DDEVertex` is a new feature describing internal delay dynamics of nodes. The resulting objects can then be delivered to the key constructor `network_dynamics`, which adds topological information of graph **`g`** and returns an `DDEFunction`, which is compatible with the solvers of `DifferentialEquations.jl`.

# Simulation
First we set the initial conditions, the time-interval `tspan`, construct the history function `h(x,p,t)` and set the parameters `p`, which are now extended to take parameters for the delay-time too.
```@example DDEVertex

x0 = randn(N) # random initial conditions
tspan = (0., 10.) #time interval

# history function defaults to all 1. and is in-place to save allocations
h(out, p, t) = (out .= 1.)

#parameters p = (vertex parameters, edge parameters, delay time)
p = (nothing, nothing, 1.)
```
Then we solve the diffusion problem.
The constructor `DDEProblem` is from the package `DealyDiffEq.jl` and wraps up the `DDEFunction` with the initial conditions `x0` , the history function `h` , time-interval `tspan` and parameters `p`, which can than be handed over to the `solve` method.

```@example DDEVertex
using OrdinaryDiffEq
using DelayDiffEq

dde_prob = DDEProblem(nd, x0, h, tspan, p)
sol = solve(dde_prob, MethodOfSteps(Tsit5()))
```
 We solve the problem with delay on the time interval $[0, 10]$ with the `Tsit5()` algorithm, recommended for solving non-stiff problems. The additional algorithm `MethodOfSteps` translates an `OrdinaryDiffEq.jl` ODE solver method into a method for delay differential equations which is highly efficient.

## Bonus: Two independet diffusions

In this extension of the first example, we will have two independent diffusions on the same network with variables $x$ and $\phi$, such that the **`dim`**=2.
First we construct the `network_dynamics`-objects.

```@example DDEVertex
nd_diffusion_vertex_2 = DDEVertex(f! = diffusionvertex!, dim = 2, sym = [:x, :ϕ])
nd_diffusion_edge_2 = StaticDelayEdge(f! = diffusionedge!, dim = 2)
nd_2 = network_dynamics(nd_diffusion_vertex_2, nd_diffusion_edge_2, g)
```
Secondly, the initial conditions are generated, where the first N values correspond to variable `x` and the values with indices from N+1 to 2N to state-variable `ϕ`, where $x \sim  N(0,1)$; $ϕ \sim N(0,1)^2$. The parameter value for the delay $\Tau$ is set to 1.0.

```@example DDEVertex
x0_2 = Array{Float64,1}(vec([randn(N).-10 randn(N).^2]')) # x ~ \mathcal{N}(0,1); ϕ ~ \mathcal{N}(0,1)^2

p = (nothing, nothing, 1.) # p = (vertexparameters, edgeparameters, delaytime)
```
Now we can define the `DDEProblem`and then solve it. 

```@Example DDEVertex
dde_prob_2 = DDEProblem(nd_2, x0_2, h, tspan, p)
sol_2 = solve(dde_prob_2, MethodOfSteps(Tsit5()));
plot(sol_2, legend=false)
```
As a solver we use again `Tsit5()` with the `MethodOfSteps`-algorithm.

## Kuramoto model with delay

An additional example which will be explained in the following is the Kuramoto model. Instead of modeling a simple diffusion on an undirected ring network, we will use the Kuramoto model for the vertices and the edges but keep delay and the Watts-Strogatz graph as in the examples above.

Again, the interactions with the neighbors are described through the edge function. Unlike the diffusion example, here the edges now have a delay. For this reason we introduce the history arrays for the destination vertices and source vertices.

```@example DDEVertex
function kuramoto_delay_edge!(e, v_s, v_d, h_v_s, h_v_d, p, t)
    # The coupling is no longer symmetric, so we need to store BOTH values (see tutorials for details)
    e[1] = p * sin(v_s[1] - h_v_d[1])
    e[2] = p * sin(h_v_s[1] - v_d[1])
    nothing
end
nothing # hide
```

Now, the contributions of the edges are summed up in each vertex so that we can define the vertex function for the Kuramoto model.

```@example DDEVertex
function kuramoto_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = p
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d
        dv[1] -= e[2]
    end
    nothing
end
```

As for the diffusion example, we now hand over the defined functions for the Kuramoto edges and vertices to the constructors `DDEVertex` and `StaticEdge` of NetworkDynamics.jl.

```@example DDEVertex
kdedge! = StaticDelayEdge(f! = kuramoto_delay_edge!, dim=2)
kdvertex! = ODEVertex(f! = kuramoto_vertex!, dim = 1)
```

Note that the edges have the dimension two, since there is no more symmetric coupling for the Kuramoto case. Accordingly, we have to set the keywork argument `dim` of `StaticDelayEdge` to two. The returned objects (`nd_diffusion_vertex` and `nd_diffusion_edge`) are passed to the key construcor `network_dynamics` together with informations of the graph which are contained in `g`.

```@example DDEVertex
nd! = network_dynamics(kdvertex!, kdedge!, g)
```

The object `nd!` returned by `network_dynamics`, is compatible with the solvers of `DifferentialEquations`.

Afterwards, random initial conditions are set for all $N$ nodes, as well as the simulation time `tspan`, the eigenfrequencies `ω` of the $N$ vertices, and the parameters `p` for the edges, the vertices, and the delay. Since we are dealing with a system containg a time delay, we need to set the history function `h` as for the examples above. The history function sets all default entries to 1.0 and is in-place to save allocations.

```@example DDEVertex
x0 = randn(N) # random initial conditions
# history function defaults to all 1. and is in-place to save allocations
h(out, p, t) = (out .= 1.)
# p = (vertexparameters, edgeparameters, delaytime)
ω = randn(N)
ω .-= sum(ω)/N
p = (ω, 2., 1.)
tspan = (0.,20.)
```

As described for the diffusion exmaple above, the constructor `DDEProblem` is used to provide an object that can be solved using the `solve` method.
