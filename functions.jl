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
