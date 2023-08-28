(* ::Package:: *)

<< FeynCalc`


iMSingleBarNu = (I Geff)/MW^2 (I g)^2 I FV[p, \[Mu]] (SpinorUBar[p3] . GA[7] . ((I GS[-q])/-SP[q]) . GA[\[Mu]] . GA[7] . SpinorV[p2, me]) (-I)/(-SP[p4 + p5] + mS^2) (SpinorU[p4] . GA[6] . SpinorV[p5])


FCClearScalarProducts[];
SP[p] = M^2; SP[p1] = Mp^2; SP[p2] = me^2; SP[p3] = SP[p4] = SP[p5] = 0;
iMSingleBarNu ComplexConjugate[iMSingleBarNu] // FermionSpinSum // DiracSimplify // FullSimplify
