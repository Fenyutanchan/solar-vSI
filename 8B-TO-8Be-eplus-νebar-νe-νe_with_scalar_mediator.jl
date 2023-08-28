# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using ProgressMeter

include("amplitude_squared_library.jl")
include("kinematics.jl")

function main(num_points::Int)

MeV = 1e6
GeV = 1e3 * MeV
mu = 931.49410242 * MeV
M = 8.0246073 * mu
Mp = 8.00530510 * mu
mS = 10 * MeV

P = [M, 0, 0, 0]
n = 5
mass_list = [Mp, 0, 0, 0, 0]

weight_sum = zeros(Real, Threads.nthreads())
prefactor = 1 / (2 * M) # for decay
prefactor /= 2 # for identical particles

p = Progress(num_points; desc="Evaluating...")
counter = Threads.Atomic{Int}(0)
ProgressMeter.update!(p, counter[])
Threads.@threads for _ ∈ 1:num_points
    ρn, pi_list = Φ(P, mass_list)
    p1, p2, p3, p4, p5 = pi_list
    q = p3 + p4 + p5

    # amplitude_squared = 4 * scalar_product(p4, p5) / (
    #     scalar_product(q)^2 *
    #     (scalar_product(p4 + p5) + mS^2)^2
    # ) * (
    #     2 * scalar_product(p3, q) * (
    #         2 * scalar_product(P, p2) * scalar_product(P, q) +
    #         M^2 * scalar_product(p2, q)
    #     ) - scalar_product(q) * (
    #         M^2 * scalar_product(p2, p3) +
    #         2 * scalar_product(P, p2) * scalar_product(P, p3)
    #     )
    # )
    amplitude_squared = amp_sq_8B_TO_8Be_e_ν̄e_νe_νe_with_scalar_mediator(P, p2, p3, p4, p5, M, mS)
    @assert amplitude_squared ≥ 0
    weight_sum[Threads.threadid()] += ρn * amplitude_squared

    Threads.atomic_add!(counter, 1)
    ProgressMeter.update!(p, counter[])
end

return sum(weight_sum) * prefactor / num_points

end

@show main(10^6)
