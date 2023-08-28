# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using CSV
using DataFrames
using LaTeXStrings
using Plots
using ProgressMeter

include("../amplitude_squared_library.jl")
include("../C_flux.jl")
include("../kinematics.jl")

function main(num_points::Int, num_plot_bar::Int)

    MeV = 1
    mu = 931.49410242 * MeV
    M = 8.0246073 * mu
    Mp = 8.00530510 * mu

    # mS_list = collect(1:2:13) * MeV
    # mS_list=[1.0, 1.5, 2, 2.5, 3. ,  4.2,  5.4,  6.6,  7.8,  9. , 10.2, 11.4, 12.6, 13.8, 15., 16, 16.5 ] * MeV
    mS_list = [10 * MeV]

    P = [M, 0, 0, 0]
    mass_list = [Mp, 0, 0, 0, 0]
    prefactor = 1 / (2 * M) # for decay

    E_ν_max = (M^2 - Mp^2) / (2 * M)
    energy_list = (collect(1:num_plot_bar) .- (1 // 2)) * E_ν_max / num_plot_bar
    flux_ν = [
        xx^2 * (M * (M - 2 * xx) - Mp^2)^2 / (M - 2 * xx)^2
        for xx ∈ energy_list
    ] * C_flux(MeV) / (8 * (2 * pi)^3)
    single_flux_νbar_collect = [zeros(num_plot_bar) for _ ∈ mS_list]
    triple_flux_νbar_collect = [zeros(num_plot_bar) for _ ∈ mS_list]
    flux_ν_max = max(flux_ν...)
    
    for (mS_index, mS) ∈ enumerate(mS_list)
        single_weighted_sum = [zeros(num_plot_bar) for _ ∈ 1:Threads.nthreads()]
        triple_weighted_sum = [zeros(num_plot_bar) for _ ∈ 1:Threads.nthreads()]

        p = Progress(num_points; desc="Evaluating @ $(mS/MeV) MeV...")
        counter = Threads.Atomic{Int}(0)
        ProgressMeter.update!(p, counter[])
        Threads.@threads for _ ∈ 1:num_points
            ρn, pi_list = Φ(P, mass_list)
            p1, p2, p3, p4, p5 = pi_list

            @assert p3[1] ≤ E_ν_max
            energy_bar_indices = [
                (Int ∘ floor)(first(p_νbar) / E_ν_max * num_plot_bar) + 1
                for p_νbar ∈ [p3, p4, p5]
            ]
    
            amplitude_squared = amp_sq_8B_TO_8Be_e_ν̄e_with_four_fermion_contact(P, p2, p3, p4, p5, M, mS)
            @assert amplitude_squared ≥ 0
            single_weighted_sum[Threads.threadid()][first(energy_bar_indices)] += ρn * amplitude_squared
            for energy_bar_index ∈ energy_bar_indices
                triple_weighted_sum[Threads.threadid()][energy_bar_index] += ρn * amplitude_squared
            end

            Threads.atomic_add!(counter, 1)
            ProgressMeter.update!(p, counter[])
        end

        x = energy_list

        ## Single
        y_νbar = sum(single_weighted_sum) * C_flux(MeV) * prefactor / num_points
        y_νbar /= (E_ν_max / num_plot_bar)

        scale_exponent = begin
            flux_νbar_max = max(y_νbar...)
            (Int ∘ floor ∘ log10)(flux_ν_max / flux_νbar_max)
        end
        y_ν = flux_ν / 10^scale_exponent

        y_min = max(
            min(filter(!iszero, y_νbar)...),
            min(filter(!iszero, y_ν)...)
        )
        y_max = max(y_νbar..., y_ν...)
        y_max = 10^((ceil ∘ log10)(y_max))

        fig = plot(
            x,
            [y_ν y_νbar];

            # title=L"Energy Spectrum of $\bar{\nu}_e$ from $^8\mathrm{B}$ Decay ($m_S$ = %$(mS / MeV) \mathrm{MeV})$",
            xlabel=L"$E_{\bar{\nu}}~[\mathrm{MeV}]$",
            # xticks=x,
            ylabel=L"Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]",
            xscale=:log10,
            yscale=:log10,
            ylims=(y_min, y_max),

            legend=:topleft,
            label=[L"$\Phi_\nu / 10^{%$scale_exponent}$" L"$\Phi_{\bar{\nu}}$"],
            minorgrid=true,
            linewidth=3
        )
        savefig(fig, "Single_νbar_Energy_Spectrum_Four-Fermion_at_$(mS/MeV)MeV.pdf")

        single_flux_νbar_collect[mS_index] = copy(y_νbar)

        ## Triple
        y_νbar = sum(triple_weighted_sum) * C_flux(MeV) * prefactor / num_points
        y_νbar /= (E_ν_max / num_plot_bar)

        scale_exponent = begin
            flux_νbar_max = max(y_νbar...)
            (Int ∘ floor ∘ log10)(flux_ν_max / flux_νbar_max)
        end
        y_ν = flux_ν / 10^scale_exponent

        y_min = max(
            min(filter(!iszero, y_νbar)...),
            min(filter(!iszero, y_ν)...)
        )
        y_max = max(y_νbar..., y_ν...)
        y_max = 10^((ceil ∘ log10)(y_max))

        fig = plot(
            x,
            [y_ν y_νbar];

            # title=L"Energy Spectrum of $\bar{\nu}_e$ from $^8\mathrm{B}$ Decay ($m_S$ = %$(mS / MeV) \mathrm{MeV})$",
            xlabel=L"$E_{\bar{\nu}}~[\mathrm{MeV}]$",
            # xticks=x,
            ylabel=L"Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]",
            xscale=:log10,
            yscale=:log10,
            ylims=(y_min, y_max),

            legend=:topleft,
            label=[L"$\Phi_\nu / 10^%$scale_exponent$" L"$\Phi_{\bar{\nu}}$"],
            minorgrid=true,
            linewidth=3
        )
        savefig(fig, "Triple_νbar_Energy_Spectrum_Four-Fermion_at_$(mS/MeV)MeV.pdf")

        triple_flux_νbar_collect[mS_index] = copy(y_νbar)
    end

    CSV.write(
        "Single_νbar_Energy_Spectrum_Four-Fermion.csv",
        DataFrame(
            "Energy"=>energy_list,
            "flux_ν"=>flux_ν,
            [
                # "flux_νbar_$(mS / MeV)MeV" => flux_νbar
                # here I (Xunjie) removed "flux_νbar_" so that it is easier to obtain mS values from the 1st line of the csv file
                "$(mS / MeV)" => flux_νbar
                for (mS, flux_νbar) ∈ zip(mS_list, single_flux_νbar_collect)
            ]...
        )
    )

    CSV.write(
        "Triple_νbar_Energy_Spectrum_Four-Fermion.csv",
        DataFrame(
            "Energy"=>energy_list,
            "flux_ν"=>flux_ν,
            [
                # "flux_νbar_$(mS / MeV)MeV" => flux_νbar
                # here I (Xunjie) removed "flux_νbar_" so that it is easier to obtain mS values from the 1st line of the csv file
                "$(mS / MeV)" => flux_νbar
                for (mS, flux_νbar) ∈ zip(mS_list, triple_flux_νbar_collect)
            ]...
        )
    )
end

main(10^7, 100)
