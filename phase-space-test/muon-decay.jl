# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using ProgressMeter

include("../kinematics.jl")

function main(num_points::Int,meMeV=0.511)

MeV = 1
GeV = 1e3 * MeV
mμ = 106 * MeV
me = meMeV * MeV
GF = 1.166e-5 * GeV^-2

M = mμ
P = [M, 0, 0, 0]
n = 3
mass_list = [0, 0, me]

weight_sum = zeros(Real, Threads.nthreads())
prefactor = 1 / (2 * M)

p = Progress(num_points; desc="Evaluating...")
counter = Threads.Atomic{Int}(0)
ProgressMeter.update!(p, counter[])
Threads.@threads for _ ∈ 1:num_points
    ρn, pi_list = Φ(P, mass_list)

    amplitude_squared = (
        64*GF^2 * scalar_product(pi_list[2], P) * scalar_product(pi_list[1], pi_list[3]) 
    )
    @assert amplitude_squared ≥ 0
    weight_sum[Threads.threadid()] += ρn * amplitude_squared

    Threads.atomic_add!(counter, 1)
    ProgressMeter.update!(p, counter[])
end

@show (GF^2 * mμ^5) / (192 * pi^3) * (1 + (8 * me^6) / mμ^6 - (8 * me^2) / mμ^2 - (12 * me^4 * log(me^2 / mμ^2)) / mμ^4)

return sum(weight_sum) * prefactor / num_points

end

@show main(10^6)

@show main(10^6, 30)
