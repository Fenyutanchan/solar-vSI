# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using LaTeXStrings
using Plots
using ProgressMeter

include("../kinematics.jl")
include("../C_flux.jl")
include("../amplitude_squared_library.jl")

function main(num_points::Int)

    MeV = 1
    mu = 931.49410242 * MeV
    M = 8.0246073 * mu
    Mp = 8.00530510 * mu

    x = collect(1:2:13)
    mS_list = x * MeV

    Γ_standard = flux_mass_factor(M, Mp) / (8 * (2 * π)^3) # * Geff^2 / MW^4
    
    P = [M, 0, 0, 0]
    n = 4
    mass_list_collect = [[Mp, 0, 0, mS] for mS ∈ mS_list]
    
    weight_sum = [
        zeros(length(mass_list_collect))
        for _ ∈ 1:Threads.nthreads()
    ]
    prefactor = 1 / (2 * M) # for decay
    
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
    
    y = sum(weight_sum) * prefactor / num_points
    y /= Γ_standard
    # y *= 100
    
    fig = plot(
        x, y;
        # title=L"Branching Ratio",
        xlabel=L"$m_S / \mathrm{MeV}$",
        xticks=x,
        ylabel=L"$R \times 100$%",
        yticks=[0.0001, 0.001, 0.01],
        # ylabel=L"$\frac{\Gamma_{^8{B} \to \mathrm{new}}}{\Gamma_{^8{B} \to \mathrm{standard}}}$ %",
        yscale=:log10,

        legend=false,
        minorgrid=true,
        linewidth=3
    )

    savefig(fig, "Branching-Ratio-Mass.pdf")

@show(x,y) # I added this line to obtain the values of x and y so that I can use them in code-xunjie/JUNO-events.ipynb
#Actually I'd prefer that (x,y) are directly written into, e.g., code-xunjie/data/branching.csv. 

end

@show main(10^6)
