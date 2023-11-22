#!/usr/bin/env wolframscript
(* ::Package:: *)

sub:={
Ct = Simplify[ct/.{Ei1 j1r1->(ej1+ej2)/2,Ei1c j1r2->(ej1-ej2)/2,e0s j0pr1->(ej3+ej4)/2,e0p j0pr2->(ej3-ej4)/2,e0sc j0pn1->(ej5+ej6)/2,e0pc j0pn2->(ej5-ej6)/2}/.
{Ei1 (j1r0-j1r3)->(ej7+ej8)/2,Ei1c (j1r0-j1r4)->(ej7-ej8)/2,Ei1 (-j1r0+j1r3)->(-ej7-ej8)/2,Ei1c (-j1r0+j1r4)->(-ej7+ej8)/2}/.
{e0s (j0pr0-j0pr3)->(ej9+ej10)/2,e0p (j0pr0-j0pr4)->(ej9-ej10)/2,e0s (-j0pr0+j0pr3)->(-ej9-ej10)/2,e0p (-j0pr0+j0pr4)->(-ej9+ej10)/2}]
}

ct = cmxA;
sub;
CmxA = Ct
ct = cmyA;
sub;
CmyA = Ct
ct = cmzA;
sub;
CmzA = Ct
ct = czxA;
sub;
CzxA = Ct
ct = czyA;
sub;
CzyA = Ct
ct = czzA;
sub;
CzzA = Ct

ct = cmxB;
sub;
CmxB = Ct
ct = cmyB;
sub;
CmyB = Ct
ct = cmzB;
sub;
CmzB = Ct
ct = czxB;
sub;
CzxB = Ct
ct = czyB;
sub;
CzyB = Ct
ct = czzB;
sub;
CzzB = Ct



