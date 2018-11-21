function R1(s,x1,x2,ϵ1)
    r1 = p.α10 + p.α11 * s + p.α12 * x1 - p.α13 * (x1^2) + p.α14 *x2 - p.α15 *(x2^2) .+ ϵ1
end

function R2(s,x1,x2,ϵ2)
    r2 = p.α20 + p.α21 * s + p.α22 * x2 - p.α23 * (x2^2) + p.α24 *x1 - p.α25 *(x1^2) .+ ϵ2
end

function R3(s,slag,ϵ3)
    if s <= 19
        r3 = p.β0 - p.β1 * (s >= 12) - p.β2 * (1 - slag) .+ ϵ3
    else
        r3 = - p.β2 .+ ϵ3 .- 40000
    end
end

function R4(ϵ4)
    r4 = p.γ0 .+ ϵ4
end

# Calculates Emax at terminal period
# SST = State Space at T

function sendto!(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, Core.eval(Main, Expr(:(=), nm, val)))
    end
end

function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto!(p; args...)
    end
end

function EMAXTi(state::Vector{Int64})
    global MC_ϵ
    cache = 0.0
    r1 = exp.(R1(state[1],state[2],state[3],MC_ϵ[1,:]))
    r2 = exp.(R2(state[1],state[2],state[3],MC_ϵ[2,:]))
    r3 = R3(state[1],state[4],MC_ϵ[3,:])
    r4 = R4(MC_ϵ[4,:])
    cache = mean(max.(r1, r2, r3, r4))
    cache
end

function EMAXTi_vec(state_vec::Vector{Vector{Int64}})
    return [EMAXTi(x) for x in state_vec]
end

function input_split(state_vec::Vector{Vector{Int64}})
    n = size(state_vec,1)
    np = nprocs()  # determine the number of processes available
    numb_splits = np-1
    x_split = Vector{Vector{eltype(state_vec)}}(undef,numb_splits)
    split_cut = convert(Int64,ceil(n/numb_splits))
    split_indx = collect(1:split_cut:n)
    for i = 1:numb_splits-1
        x_split[i] = state_vec[split_indx[i]:(split_indx[i+1]-1)]
    end
    x_split[numb_splits] = state_vec[split_indx[numb_splits]:n]
    return x_split
end

function output_vcat(output_vcat)
    expr_str = ""
    expr_str = expr_str*"vcat(Emax[1]"
    for i=2:numb_splits
        expr_str = expr_str * ",Emax[$i]"
    end
    expr_str =  expr_str*")"
    return Meta.parse(expr_str)
end

function pmap_EMAXTi(state_vec::Vector{Vector{Int64}})
    np = nprocs()  # determine the number of processes available
    x_split = input_split(state_vec)
    n = size(x_split,1)
    numb_splits = n
    Emax = Array{Array{Float64}}(undef,numb_splits)
    i = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for processor=1:np
            if processor != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
                        Emax[idx] = remotecall_fetch(EMAXTi_vec,processor,x_split[idx])
                    end
                end
            end
        end
    end
    vcat(Emax[1],Emax[2],Emax[3],Emax[4],Emax[5],Emax[6],Emax[7])
end

function EMAXti(state::Vector{Int64})
    global fEmax
    global MC_ϵ
    d3 = state + [1,0,0,1-state[4]]
    d3[1] = min(20,d3[1])
    v1 = exp.(R1(state[1],state[2],state[3],MC_ϵ[1,:])) .+ p.β*fEmax[state+[0,1,0,-state[4]]]
    v2 = exp.(R2(state[1],state[2],state[3],MC_ϵ[2,:])) .+ p.β*fEmax[state+[0,0,1,-state[4]]]
    if d3[1] < 20
        v3 = R3(state[1],state[4],MC_ϵ[3,:]) .+ p.β*fEmax[d3]
    else
        v3 = R3(state[1],state[4],MC_ϵ[3,:])
    end
    v4 = R4(MC_ϵ[4,:]) .+ p.β*fEmax[state+[0,0,0,-state[4]]]
    cache = mean(max.(v1, v2, v3, v4))
    cache
end

function EMAXti_vec(state_vec::Vector{Vector{Int64}})
    return [EMAXti(x) for x in state_vec]
end

function pmap_EMAXti(state_vec::Vector{Vector{Int64}})
    n = size(state_vec,1)
    np = nprocs()  # determine the number of processes available
    if np==1 || n<100
        Emax = remotecall_fetch(EMAXti_vec,1,state_vec)
        return Emax
    else
        x_split = input_split(state_vec)
        n = size(x_split,1)
        numb_splits = n
        Emax = Array{Array{Float64}}(undef,numb_splits)
        i = 1
        # function to produce the next work item from the queue.
        # in this case it's just an index.
        nextidx() = (idx=i; i+=1; idx)
        @sync begin
            for processor=1:np
                if processor != myid() || np == 1
                    @async begin
                        while true
                            idx = nextidx()
                            if idx >n
                                break
                            end
                            Emax[idx] = remotecall_fetch(EMAXti_vec,processor,x_split[idx])
                        end
                    end
                end
            end
        end
        return vcat(Emax[1],Emax[2],Emax[3],Emax[4],Emax[5],Emax[6],Emax[7])
    end
end
