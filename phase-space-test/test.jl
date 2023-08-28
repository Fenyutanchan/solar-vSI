# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Dates
using ProgressMeter
using Random
using SHA

set_random_seed_from_time() = 
    (Random.seed! ∘ Int ∘ mod)(
        parse(BigInt,
            (bytes2hex ∘ sha256 ∘ string ∘ now)(),
            base=16
        ),
        typemax(Int)
    )

function main(num_points::Int)
    p = Progress(num_points; desc="Generating random points")
    counter = Threads.Atomic{Int}(0)
    ProgressMeter.update!(p, counter[])

    inv_num_points = inv(num_points)

    result1 = zeros(Real, Threads.nthreads())
    result2 = zeros(Real, Threads.nthreads())
    Threads.@threads for _ ∈ 1:num_points
        set_random_seed_from_time()
        x1, x2 = rand(2)
        x1 + x2 < 1 && (result1[Threads.threadid()] += inv_num_points)

        x1 = rand()
        x2 = rand() * (1 - x1)
        result2[Threads.threadid()] += (1 - x1) * inv_num_points

        Threads.atomic_add!(counter, 1)
        ProgressMeter.update!(p, counter[])
    end

    return sum(result1), sum(result2)
end

@show main(10^7)
