# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Dates
using Random
using SHA

λ(x, y, z) = (x - y - z)^2 - 4 * y * z
set_random_seed_from_time() = 
    (Random.seed! ∘ Int ∘ mod)(
        parse(BigInt,
            (bytes2hex ∘ sha256 ∘ string ∘ now)(),
            base=16
        ),
        typemax(Int)
    )
function scalar_product(p::Vector{<:Real}, q::Vector{<:Real})::Real
    @assert length(p) == length(q) == 4
    return p[1] * q[1] - p[2] * q[2] - p[3] * q[3] - p[4] * q[4]
end
scalar_product(p) = scalar_product(p, p)
calc_mass(p) = (sqrt∘scalar_product)(p)

function calc_lorentz_boost_matrix(p::Vector{<:Real})::Matrix{Real}
    @assert length(p) == 4
    m = calc_mass(p)
    @assert !iszero(m) "Invalid: massless four-momentum $p."
    γ = p[1] / m
    pp = p[2]^2 + p[3]^2 + p[4]^2
    if iszero(pp)
        @assert iszero(γ - 1)
        return [
            1 0 0 0
            0 1 0 0
            0 0 1 0
            0 0 0 1
        ]
    end

    Λ = Matrix{Real}(undef, 4, 4)
    Λ[1, 1] = γ
    for ii ∈ 2:4
        Λ[1, ii] = Λ[ii, 1] = p[ii] * sqrt((γ^2 - 1) / pp) 
    end
    for ii ∈ 2:4, jj ∈ ii:4
        Λ[ii, jj] = Λ[jj, ii] = (ii == jj ? 1 : 0) + (γ - 1) * p[ii] * p[jj] / pp
    end
    return Λ
end
function lorentz_boost(p::Vector{<:Real}, ref_momentum::Vector{<:Real})::Vector{Real}
    @assert length(p) == length(ref_momentum) == 4
    Λ = calc_lorentz_boost_matrix(ref_momentum)
    return Λ * p
end

function Φ(
    P::Vector{<:Real},
    # n::Int,
    mass_list::Vector{<:Real}
)::Tuple{Real, Vector{Vector{<:Real}}}
    @assert length(P) == 4 "Invalid four-momentum $P."
    # @assert length(mass_list) == n "The number of masses is not equal to $n."
    n = length(mass_list)

    set_random_seed_from_time()

    M_list = Vector{Real}(undef, n)
    M_list[1] = mass_list[1]
    M_list[n] = calc_mass(P)
    ρn = 1
    # while true
    #     M_condition = true
    #     for ii ∈ n-1:-1:2
    #         μi = sum(mass_list[1:ii])
    #         integration_range = M_list[ii+1] - sum(mass_list[ii+1:end]) - μi
    #         M_list[ii] = rand() * integration_range .+ μi
    #         ρn *= integration_range
    #         if M_list[ii] ≥ M_list[ii+1] - mass_list[ii+1]
    #             M_condition = false
    #         end
    #     end
    #     M_condition && break
    # end
    for ii ∈ n-1:-1:2
        μi = sum(mass_list[1:ii])
        integration_range = M_list[ii+1] - mass_list[ii+1] - μi
        M_list[ii] = rand() * integration_range .+ μi
        ρn *= integration_range
    end
    ρ_list = [
        (sqrt ∘ λ)(
            M_list[ii]^2,
            M_list[ii-1]^2,
            mass_list[ii]^2
        ) / (2 * M_list[ii])
        for ii ∈ 2:n
    ]
    ρn *= prod(ρ_list) / (2^n * M_list[n] * (2 * pi)^(3 * n - 4))

    pi_list = Vector{Vector{Real}}(undef, n)
    ki_list = Vector{Vector{Real}}(undef, n)
    ki_list[n] = P
    θ_list = acos.(rand(n) * 2 .- 1)
    ϕ_list = rand(n) * 2 * π
    ρn *= (4 * pi)^(n - 1)
    
    for ii ∈ n:-1:2
        Ei = (M_list[ii]^2 + mass_list[ii]^2 - M_list[ii-1]^2) / (2 * M_list[ii])
        pp = sqrt(Ei^2 - mass_list[ii]^2)
        Ek = (M_list[ii]^2 + M_list[ii-1]^2 - mass_list[ii]^2) / (2 * M_list[ii])

        rest_pi = [
            Ei,
            pp * sin(θ_list[ii]) * cos(ϕ_list[ii]),
            pp * sin(θ_list[ii]) * sin(ϕ_list[ii]),
            pp * cos(θ_list[ii])
        ]
        pi_list[ii] = lorentz_boost(rest_pi, ki_list[ii])

        rest_ki = [Ek, -rest_pi[2:4]...]
        ki_list[ii-1] = lorentz_boost(rest_ki, ki_list[ii])
    end
    pi_list[1] = ki_list[1]

    return ρn, pi_list
end

function test_kinematics()
    MeV = 1
    mu = 931.49410242 * MeV
    M = 8.0246073 * mu
    Mp = 8.00530510 * mu

    P = [M, 0, 0, 0]
    n = 3
    mass_list = [Mp, 0, 0]

    ρn, pi_list = Φ(P, mass_list)

    @assert sum(pi_list) ≈ P
    
    return ρn, pi_list
end

# test_kinematics()
