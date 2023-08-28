# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using ProgressMeter

include("../C_flux.jl")
include("../amplitude_squared_library.jl")
include("../kinematics.jl")

function main(num_points::Int)

MeV = 1
# GeV = 1e3 * MeV
mu = 931.49410242 * MeV
M = 8.0246073 * mu
Mp = 8.00530510 * mu
# Mp = 1.00727647 * mu

P = [M, 0, 0, 0]
n = 3
mass_list = [Mp, 0, 0]

weight_sum = zeros(Real, Threads.nthreads())
prefactor = 1 / (2 * M)

p = Progress(num_points; desc="Evaluating...")
counter = Threads.Atomic{Int}(0)
ProgressMeter.update!(p, counter[])
Threads.@threads for _ ∈ 1:num_points
    ρn, pi_list = Φ(P, mass_list)

    # amplitude_squared = (
    #     4 * scalar_product(pi_list[2], P) * scalar_product(pi_list[3], P) -
    #     2 * scalar_product(pi_list[2], pi_list[3]) * scalar_product(P)
    # )
    amplitude_squared = amp_sq_8B_TO_8Be_e_νe(P, pi_list[2], pi_list[3])
    @assert amplitude_squared ≥ 0
    weight_sum[Threads.threadid()] += ρn * amplitude_squared

    Threads.atomic_add!(counter, 1)
    ProgressMeter.update!(p, counter[])
end

# Wrong expression
# @show (
#     (M^2 - Mp^2) * (M^4 + 10 * M^2 * Mp^2 + Mp^4) / (24 * M) +
#     M * Mp^2 * (M^2 + Mp^2) * log(Mp / M) / 2
# ) / (8 * (2 * pi)^3)

# Correct expression
# @show (
#     M^8 - 8 * M^6 * Mp^2 + 24 * M^4 * Mp^4 * log(M / Mp)
#     + 8 * M^2 * Mp^6 - Mp^8
# ) / (96 * M^3) / (8 * (2 * pi)^3)
@show flux_mass_factor(M, Mp) / (8 * (2 * pi)^3)

return sum(weight_sum) * prefactor / num_points

end

@show main(10^7)
