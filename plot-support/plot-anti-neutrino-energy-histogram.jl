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

    mS_list = collect(1:2:13) * MeV
    mS_list=[1.0, 1.5, 2, 2.5, 3. ,  4.2,  5.4,  6.6,  7.8,  9. , 10.2, 11.4, 12.6, 13.8, 15., 16, 16.5 ] * MeV
    # mS_list = [6 * MeV]

    P = [M, 0, 0, 0]
    n = 4
    mass_list_collect = [[Mp, 0, 0, mS] for mS ∈ mS_list]
    prefactor = 1 / (2 * M) # for decay

    E_ν_max = (M^2 - Mp^2) / (2 * M)
    energy_list = (collect(1:num_plot_bar) .- (1 // 2)) * E_ν_max / num_plot_bar
    flux_ν = [
        xx^2 * (M * (M - 2 * xx) - Mp^2)^2 / (M - 2 * xx)^2
        for xx ∈ energy_list
    ] * C_flux(MeV) / (8 * (2 * pi)^3)
    flux_νbar_collect = [zeros(num_plot_bar) for _ ∈ mS_list]
    flux_ν_max = max(flux_ν...)
    
    for (mass_index, mass_list) ∈ enumerate(mass_list_collect)
        weighted_sum = [zeros(num_plot_bar) for _ ∈ 1:Threads.nthreads()]
        additional_weighted_sum = [zeros(num_plot_bar) for _ ∈ 1:Threads.nthreads()]
        mS = last(mass_list)
        E_νbar_max = begin
            original_E_νbar_max = (M^2 - (Mp + mS)^2) / (2 * M)
            E_Φ_max = (M^2 - Mp^2) / (2 * M)
    
            γ = E_Φ_max / mS
            additional_E_νbar_max = (sqrt(γ^2 - 1) + γ) * mS / 2

            max(original_E_νbar_max, additional_E_νbar_max)
        end     

        p = Progress(num_points; desc="Evaluating @ $mS MeV...")
        counter = Threads.Atomic{Int}(0)
        ProgressMeter.update!(p, counter[])
        Threads.@threads for _ ∈ 1:num_points
            ρn, pi_list = Φ(P, mass_list)
            p1, p2, p3, p4 = pi_list

            @assert p3[1] ≤ E_νbar_max
            energy_bar_index = (Int ∘ floor)(p3[1] / E_ν_max * num_plot_bar) + 1
    
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
            weighted_sum[Threads.threadid()][energy_bar_index] += ρn * amplitude_squared

            _, additional_pi_list = Φ(p4, [0, 0])
            for additional_pi ∈ additional_pi_list
                energy_bar_index = (Int ∘ floor)(additional_pi[1] / E_ν_max * num_plot_bar) + 1
                additional_weighted_sum[Threads.threadid()][energy_bar_index] += ρn * amplitude_squared / 2
            end

            Threads.atomic_add!(counter, 1)
            ProgressMeter.update!(p, counter[])
        end

        x = energy_list
        original_y_νbar = sum(weighted_sum) * C_flux(MeV) * prefactor / num_points
        final_y_νbar = original_y_νbar + sum(additional_weighted_sum) * C_flux(MeV) * prefactor / num_points
        original_y_νbar /= (E_ν_max / num_plot_bar)
        final_y_νbar /= (E_ν_max / num_plot_bar)

        scale_exponent = begin
            flux_νbar_max = max(original_y_νbar..., final_y_νbar...)
            (Int ∘ floor ∘ log10)(flux_ν_max / flux_νbar_max)
        end
        y_ν = flux_ν / 10^scale_exponent

        y_min = max(
            min(filter(!iszero, original_y_νbar)...),
            min(filter(!iszero, final_y_νbar)...),
            min(filter(!iszero, y_ν)...)
        )
        y_max = max(original_y_νbar..., final_y_νbar..., y_ν...)
        y_max = 10^((ceil ∘ log10)(y_max))

        fig = plot(
            x,
            [y_ν original_y_νbar final_y_νbar];

            # title=L"Energy Spectrum of $\bar{\nu}_e$ from $^8\mathrm{B}$ Decay ($m_S$ = %$(mS / MeV) \mathrm{MeV})$",
            xlabel=L"$E_{\bar{\nu}}~[\mathrm{MeV}]$",
            # xticks=x,
            ylabel=L"Flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]",
            xscale=:log10,
            yscale=:log10,
            ylims=(y_min, y_max),

            legend=:topleft,
            label=[L"$\Phi_\nu / 10^{%$scale_exponent}$" L"$\Phi_{\bar{\nu}}$ (original)" L"$\Phi_{\bar{\nu}}$ (final)"],
            minorgrid=true,
            linewidth=3
        )
        savefig(fig, "Energy_Spectrum_at_$(mS / MeV)MeV.pdf")

        flux_νbar_collect[mass_index] = copy(final_y_νbar)
        # while length(flux_νbar) < length(energy_list)
        #     push!(flux_νbar, 0)
        # end
    end

    CSV.write(
        "Energy_Spectrum.csv",
        DataFrame(
            "Energy"=>energy_list,
            "flux_ν"=>flux_ν,
            [
                # "flux_νbar_$(mS / MeV)MeV" => flux_νbar
                # here I (Xunjie) removed "flux_νbar_" so that it is easier to obtain mS values from the 1st line of the csv file
                "$(mS / MeV)" => flux_νbar
                for (mS, flux_νbar) ∈ zip(mS_list, flux_νbar_collect)
            ]...
        )
    )
end

# main(10^7, 100)
# it turns out 10^5 is enough to guarantee that the final result is a smooth curve.
main(10^5, 100)
