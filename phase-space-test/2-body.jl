# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using ProgressMeter

include("../kinematics.jl")

function main(num_points::Int)

P = [10, 0, 0, 0]
n = 2
mass_list = [0, 0]

weight_sum = zeros(Real, Threads.nthreads())
prefactor = 1

p = Progress(num_points; desc="Evaluating...")
counter = Threads.Atomic{Int}(0)
ProgressMeter.update!(p, counter[])
Threads.@threads for _ ∈ 1:num_points
    ρn, pi_list = Φ(P, mass_list)

    weight_sum[Threads.threadid()] += ρn

    Threads.atomic_add!(counter, 1)
    ProgressMeter.update!(p, counter[])
end

@show 5 / (4 * pi * 10)

return sum(weight_sum) * prefactor / num_points

end

@show main(10^7)
