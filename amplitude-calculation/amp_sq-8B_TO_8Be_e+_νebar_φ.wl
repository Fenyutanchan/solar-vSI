(* ::Package:: *)

<< FeynCalc`


FCClearScalarProducts[];
SP[p] = M^2; SP[p2] = me^2; SP[p3] = 0; SP[p4] = mS^2;
iM = (I Geff)/MW^2 I g I FV[p, \[Mu]] (SpinorUBar[p3] . GA[7] . ((I GS[p3 + p4])/SP[p3 + p4]) . GA[\[Mu]] . GA[7] . SpinorV[p2, me])
ComplexConjugate[iM] iM //FermionSpinSum //DiracSimplify //ScalarProductExpand //FullSimplify



