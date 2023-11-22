#!/usr/bin/env wolframscript
(* ::Package:: *)

sub:={
Ct = Simplify[ct/.{Eis jsr1->(ej1+ej2)/2,Eisc jsr2->(ej1-ej2)/2,e0s j0nr1->(ej3+ej4)/2,e0p j0nr2->(ej3-ej4)/2}/.
{Eis (jsr0-jsr3)->(ej5+ej6)/2,Eisc (jsr0-jsr4)->(ej5-ej6)/2,Eis (-jsr0+jsr3)->(-ej5-ej6)/2,Eisc (-jsr0+jsr4)->(-ej5+ej6)/2}/.
{e0s (j0nr0-j0nr3)->(ej7+ej8)/2,e0p (j0nr0-j0nr4)->(ej7-ej8)/2,e0s (-j0nr0+j0nr3)->(-ej7-ej8)/2,e0p (-j0nr0+j0nr4)->(-ej7+ej8)/2}]
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



