#!/usr/bin/env wolframscript
(* ::Package:: *)

sub:={
Ct = Simplify[ct/.{Eis jsr1->(ej1+ej2)/2,Eisc jsr2->(ej1-ej2)/2,e0sc j0rn1->(ej3+ej4)/2,e0pc j0rn2->(ej3-ej4)/2}/.
{Eis (jsr0-jsr3)->(ej5+ej6)/2,Eisc (jsr0-jsr4)->(ej5-ej6)/2,Eis (-jsr0+jsr3)->(-ej5-ej6)/2,Eisc (-jsr0+jsr4)->(-ej5+ej6)/2}/.
{e0sc (j0rn0-j0rn3)->(ej7+ej8)/2,e0pc (j0rn0-j0rn4)->(ej7-ej8)/2,e0sc (-j0rn0+j0rn3)->(-ej7-ej8)/2,e0pc (-j0rn0+j0rn4)->(-ej7+ej8)/2}]
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

