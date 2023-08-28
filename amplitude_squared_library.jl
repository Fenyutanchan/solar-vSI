# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

include("kinematics.jl")

# ⁸B → ⁸Be + e⁺ + νₑ
amp_sq_8B_TO_8Be_e_νe(P, p2, p3) = (
    4 * scalar_product(p2, P) * scalar_product(p3, P) -
    2 * scalar_product(p2, p3) * scalar_product(P)
) # * Geff^2 / MW^4

# ⁸B → ⁸Be + e⁺ + ν̄ₑ + φ
# Note: the expression used in amp_sq_8B_TO_8Be_e_ν̄e_φ is equivalent to this one:
# (scalar_product(p,p2)*(-4*q2*scalar_product(p,p3) + 8*scalar_product(p,q)*scalar_product(p3,q)) + 2*M^2*(q2*scalar_product(p2,p3) - 2*scalar_product(p2,q)*scalar_product(p3,q)))/q2^2
amp_sq_8B_TO_8Be_e_ν̄e_φ(P, p2, p3, p4, M, mS) = -2 / (
    2 * scalar_product(p3, p4) + mS^2
)^2 * (
    2 * M^2 * scalar_product(p2, p4) * scalar_product(p3, p4)
    - M^2 * mS^2 * scalar_product(p2, p3)
    - 4 * scalar_product(p2, P) * scalar_product(p4, P) * scalar_product(p3, p4)
    + 2 * mS^2 * scalar_product(p2, P) * scalar_product(p3, P)
) # * g^2 * Geff^2 / MW^4
    
function amp_sq_8B_TO_8Be_e_ν̄e_νe_νe_with_scalar_mediator(P, p2, p3, p4, p5, M, mS)
    q = p3 + p4 + p5
    return 4 * scalar_product(p4, p5) / (
        scalar_product(q)^2 * (
            scalar_product(p4 + p5) - mS^2
        )^2
    ) * (
        2 * scalar_product(p3, q) * (
            2 * scalar_product(p2, P) * scalar_product(q, P)
            - M^2 * scalar_product(p2, q)
        ) - scalar_product(q) * (
            -M^2 * scalar_product(p2, p3)
            + 2 * scalar_product(p2, P) * scalar_product(p3, P)
        )
    ) # * g^2 * Geff^2 / MW^4
end

# ⁸B → ⁸Be + e⁺ + ν̄ₑ + ( νₑ + νₑ or ν̄ₑ + ν̄ₑ ), four-fermion contact interaction
function amp_sq_8B_TO_8Be_e_ν̄e_with_four_fermion_contact(P, p2, p3, p4, p5, M, mS)
    q = p3 + p4 + p5
    return 4 * scalar_product(p4, p5) / (
        scalar_product(q)^2 * mS^4
    ) * (
        2 * scalar_product(p3, q) * (
            2 * scalar_product(p2, P) * scalar_product(q, P)
            - M^2 * scalar_product(p2, q)
        ) - scalar_product(q) * (
            -M^2 * scalar_product(p2, p3)
            + 2 * scalar_product(p2, P) * scalar_product(p3, P)
        )
    ) # * Geff_4f^2 * Geff^2 / MW^4
end

function amp_sq_S(p, p2, p3, p4, p5, M)
    q = p3 + p4 + p5
    q2=scalar_product(q)

    # expression automatically generated via xxx.nb
    return (-4*(2*scalar_product(p,p2)*(q2*scalar_product(p,p3) - 2*scalar_product(p,q)*scalar_product(p3,q)) + M^2*(-(q2*scalar_product(p2,p3)) + 2*scalar_product(p2,q)*scalar_product(p3,q)))*scalar_product(p4,p5))/q2^2 # * g_eff^2 / MW^4 / Λ^4
end


function amp_sq_V(p, p2, p3, p4, p5, M)
    q = p3 + p4 + p5
    q2=scalar_product(q)

    # expression automatically generated via xxx.nb
    return (-16*scalar_product(p3,p4)*(2*scalar_product(p,p2)*(q2*scalar_product(p,p5) - 2*scalar_product(p,q)*scalar_product(p5,q)) + M^2*(-(q2*scalar_product(p2,p5)) + 2*scalar_product(p2,q)*scalar_product(p5,q))))/q2^2 # * g_eff^2 / MW^4 / Λ^4
end


function amp_sq_T(p, p2, p3, p4, p5, M)
    q = p3 + p4 + p5
    q2=scalar_product(q)

    # expression automatically generated via xxx.nb
    return (-64*(scalar_product(p,p2)*(4*q2*scalar_product(p,p5)*scalar_product(p3,p4) + 4*q2*scalar_product(p,p4)*scalar_product(p3,p5) - 2*q2*scalar_product(p,p3)*scalar_product(p4,p5) + 4*scalar_product(p,q)*scalar_product(p3,q)*scalar_product(p4,p5) - 8*scalar_product(p,q)*scalar_product(p3,p5)*scalar_product(p4,q) - 8*scalar_product(p,q)*scalar_product(p3,p4)*scalar_product(p5,q)) + M^2*(-2*q2*scalar_product(p2,p5)*scalar_product(p3,p4) - 2*q2*scalar_product(p2,p4)*scalar_product(p3,p5) + q2*scalar_product(p2,p3)*scalar_product(p4,p5) - 2*scalar_product(p2,q)*scalar_product(p3,q)*scalar_product(p4,p5) + 4*scalar_product(p2,q)*scalar_product(p3,p5)*scalar_product(p4,q) + 4*scalar_product(p2,q)*scalar_product(p3,p4)*scalar_product(p5,q))))/q2^2 # * g_eff^2 / MW^4 / Λ^4
end


function amp_sq_ϕ(p, p2, p3, p4, M)
    q = p3 + p4 
    q2=scalar_product(q)

    # expression automatically generated via xxx.nb
    return (scalar_product(p,p2)*(-4*q2*scalar_product(p,p3) + 8*scalar_product(p,q)*scalar_product(p3,q)) + 2*M^2*(q2*scalar_product(p2,p3) - 2*scalar_product(p2,q)*scalar_product(p3,q)))/q2^2 # * g_eff^2 / MW^4 
end