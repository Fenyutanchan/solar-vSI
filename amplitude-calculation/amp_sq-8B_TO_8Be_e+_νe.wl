(* ::Package:: *)

<< FeynCalc`


iM = (I Geff)/MW^2 I FV[p, \[Mu]] (SpinorUBar[k3] . GA[\[Mu]] . GA[7] . SpinorV[k2]);
FCClearScalarProducts[];

Msqr = ComplexConjugate[iM] iM //FermionSpinSum //DiracSimplify //FullSimplify
SP[p] = M^2; SP[k2, p] = M E2; SP[k3, p] = M E3; SP[k2, k3] = E2 E3 (1 - Cos[\[Theta]23]);
Msqr


q = ((M-E2-E3)^2 - Mp^2 - E2^2 - E3^2) / (2 E2 E3) //FullSimplify;
1 + q //FullSimplify //Numerator
FullSimplify@Integrate[%, {E2, (M (M - 2 E3) - Mp^2) / (2M), (M (M - 2 E3) - Mp^2) / (2 (M - 2E3))}]


Solve[(M - E3 - E2)^2 == Mp^2 + E3^2 + E2^2 + 2 E2 E3 q23, E2]


LogLogPlot[-((E3^2 (2 E3 M - M^2 + Mp^2)^2)/(2 E3 - M)),{E3, 0, (M^2 - Mp^2) / (2 M)}] /.{M -> 10, Mp -> 9}


int = Integrate[-((E3^2 (2 E3 M-M^2+Mp^2)^2)/(M(2 E3-M))), {E3, 0, (M^2 - Mp^2) / (2 M)}, Assumptions -> M > Mp > 0] //FullSimplify[#, Assumptions -> M > Mp > 0]&


Cflux = \[CapitalPhi] / (int / (4 \[Pi] AU^2 64 \[Pi]^3)) //.{M -> 8.0246073 mu, Mp -> 8.00530510 mu, mu -> 931.49410242 MeV, MeV -> 10^-3 GeV, AU -> 1.495978707*^11 m, \[CapitalPhi] -> 5.16*^6 cm^-2 s^-1, cm -> 10^-2 m, s->(6.582119569*^-22 MeV)^-1}
