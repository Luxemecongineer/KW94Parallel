# Orginal script is coded by Marcos Lee.
using Distributed
using Random
using Distributions
using LinearAlgebra

@everywhere import DataStructures.OrderedDict
@everywhere import Statistics: mean


@everywhere include("feasibleSet.jl")
@everywhere include("functions.jl")
@everywhere param = 1       #which parameter set to use
@everywhere include("Parameters$param.jl")

include("function_1stCore.jl")

N = 1000        #Number of people
T = 40          #I start with t =1, the paper starts with t = 0
MC = 2000     #Number of MC draws
st = [10,0,0,1] # Initial state at t=1

Domain_set = StateSpace(st, T)

mu = [0, 0, 0, 0] #Mean of ϵ
sigma = [p.σ11^2 p.σ12 p.σ13 p.σ14;
        p.σ12 p.σ22^2 p.σ23 p.σ24;
        p.σ13 p.σ23 p.σ33^2 p.σ34;
        p.σ14 p.σ24 p.σ34 p.σ44^2]


Random.seed!(10)
MC_ϵ = rand(MvNormal(mu, sigma),MC) #Take S draws from the Multivariate Normal
@time Emaxall = genEmaxAll(Domain_set,MC_ϵ, T)
