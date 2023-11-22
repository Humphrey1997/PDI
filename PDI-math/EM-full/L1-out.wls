#!/usr/bin/env wolframscript
(* ::Package:: *)

expand :={
terms = {};
count = 0;
Ta = T0;
indexa = StringPart[ToString[a],2;;-1];
indexpa = "1";
cu = Simplify[CoefficientList[Ta,u]];
numu = Part[Dimensions[cu],1];
Do[
term = Simplify[Part[cu,iu]];
If[ToString[term]!="0",
indexpu = ToString[iu-1];
Zg = fg<>indexpu<>indexpa;
Zf = ff<>ToString[count+1];
count = count + 1;
G = "g1("<>indexpu<>",zeta"<>indexa<>","<>indexpa<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa;
Zff = Zff<>Zf<>"="<>Zg<>";";
Zgg = Zgg<>Zg<>"="<>G<>CS<>";\n";
terms = Append[terms,term];
,Print["Jump"]];
,{iu,numu}];

Zff = Zff<>"\n";
}

a=b1n;pa=1;
DC = a^pa;
Zff = {};
Zgg = {};

fg = "ZA";

T0 = cmxA*DC;
ff = "ZAmx";
expand;
termsAmx = terms;

T0 = cmyA*DC;
ff = "ZAmy";
expand;
termsAmy = terms;

T0 = cmzA*DC;
ff = "ZAmz";
expand;
termsAmz = terms;

T0 = czxA*DC;
ff = "ZAzx";
expand;
termsAzx = terms;

T0 = czyA*DC;
ff = "ZAzy";
expand;
termsAzy = terms;

T0 = czzA*DC;
ff = "ZAzz";
expand;
termsAzz = terms;

Print[Zgg]
Print[Zff]



