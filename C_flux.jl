# Copyright (c) 2023 Quan-feng WU
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

flux_mass_factor(M, Mp) = (
    M^8 - 8 * M^6 * Mp^2 + 24 * M^4 * Mp^4 * log(M / Mp)
    + 8 * M^2 * Mp^6 - Mp^8
) / (96 * M^3)

function C_flux(MeV)
    mu = 931.49410242 * MeV
    M = 8.0246073 * mu
    Mp = 8.00530510 * mu
    # s = (6.582119569e-22 * MeV)^-1

    constant_factor = 1 / (64 * π^3)
    
    flux = 5.16e6 # * cm^-2 * s^-1

    return flux / (constant_factor * flux_mass_factor(M, Mp))
end

function test_C_flux()
    MeV = 10^-3
    m = 10^2
    AU = 1.495978707e11 * m
    s = (6.582119569e-22 * MeV)^-1

    return C_flux(MeV) * (4 * π * AU^2) * s^-1 # in GeV^-4
end

# @show test_C_flux()
