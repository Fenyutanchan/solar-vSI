# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using ProgressMeter

include("amplitude_squared_library.jl")
include("C_flux.jl")
include("kinematics.jl")

function main(num_points::Int)

cm = 1
MeV = 1
GeV = 1e3 * MeV
mu = 931.49410242 * MeV
M = 8.0246073 * mu
Mp = 8.00530510 * mu
mS_list = collect(1:2:13) * MeV
# s = (6.582119569e-22 * MeV)^-1
factor_MeV_to_s_inv = 6.582119569e-22^-1

P = [M, 0, 0, 0]
n = 4
mass_list_collect = [[Mp, 0, 0, mS] for mS ∈ mS_list]

weight_sum = [
    zeros(length(mass_list_collect))
    for _ ∈ 1:Threads.nthreads()
]
prefactor = 1 / (2 * M) # for decay
# prefactor /= 2 # for identical particles

p = Progress(num_points; desc="Evaluating...")
counter = Threads.Atomic{Int}(0)
ProgressMeter.update!(p, counter[])
Threads.@threads for _ ∈ 1:num_points
    for (mass_index, mass_list) ∈ enumerate(mass_list_collect)
        ρn, pi_list = Φ(P, mass_list)
        mS = last(mass_list)
        p1, p2, p3, p4 = pi_list

        # Wrong expression
        # amplitude_squared = 1 / (
        #     2 * (2 * scalar_product(p3, p4) - mS^2)^2
        # ) * (
        #     2 * M^2 * scalar_product(p2, p4) * scalar_product(p3, p4)
        #     - M^2 * mS^2 * scalar_product(p2, p3)
        #     - 4 * scalar_product(P, p2) * scalar_product(P, p4) * scalar_product(p3, p4)
        #     + 2 * mS^2 * scalar_product(P, p2) * scalar_product(P, p3)
        # )
        amplitude_squared = amp_sq_8B_TO_8Be_e_ν̄e_φ(P, p2, p3, p4, M, mS)
        @assert amplitude_squared ≥ 0
        weight_sum[Threads.threadid()][mass_index] += ρn * amplitude_squared
    end

    Threads.atomic_add!(counter, 1)
    ProgressMeter.update!(p, counter[])
end

return C_flux(MeV) * sum(weight_sum) * prefactor / num_points
# return sum(weight_sum) * prefactor / num_points

end

@show main(10^8)
