#!/usr/bin/env wolframscript
(* ::Package:: *)

a=b0n;pa=1;
DC = a^pa;
N0 = NC*DC;
Zff = "";
count = 0;
termsA = {};
ca = CoefficientList[N0,a];
numa = Part[Dimensions[ca],1];
Do[
TA = Part[ca,ia];
cu = Simplify[CoefficientList[TA,u]];
numu = Part[Dimensions[cu],1];
Do[
term = Part[cu,iu];
If[ToString[term]!="0",
indexa = StringPart[ToString[a],2;;-1];
indexpa = ToString[pa-ia+1];
indexu = ToString[iu-1];
Zf = "ZA"<>ToString[count+1];
count = count + 1;
termsA = Append[termsA,term];
G = "g1("<>indexu<>",zeta"<>indexa<>","<>indexpa<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa;
Zff = Zff<>Zf<>"="<>G<>CS<>";\n";
,Print["Jump"]];
,{iu,numu}];
,{ia,numa}];
Print[Zff]
