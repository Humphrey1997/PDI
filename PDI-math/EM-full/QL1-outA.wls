#!/usr/bin/env wolframscript
(* ::Package:: *)

expand :={
terms = {};
count = 0;
Ta = T0;
indexa = StringPart[ToString[a],2;;-1];
indexpa = "1";
cb = CoefficientList[Ta,b];
numb = Part[Dimensions[cb],1];
Do[
Tb = Part[cb,ib];
cu = Simplify[CoefficientList[Tb,u]];
numu = Part[Dimensions[cu],1];
Do[
term = FullSimplify[Part[cu,iu]];
If[ToString[term]!="0",
indexb = StringPart[ToString[b],2;;-1];
indexpb = ToString[pb-ib+1];
indexpu = ToString[iu-1];
Zg = fg<>indexpu<>indexpa<>indexpb;
Zf = ff<>ToString[count+1];
count = count + 1;
G = "g2("<>indexpu<>",zeta"<>indexa<>","<>indexpa<>",zeta"<>indexb<>","<>indexpb<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa<>"/c"<>StringPart[ToString[b],2]<>"^"<>indexpb;
Zff = Zff<>Zf<>"="<>Zg<>";";
Zgg = Zgg<>Zg<>"="<>G<>CS<>";\n";
terms = Append[terms,term];
,Print["Jump"]];
,{iu,numu}];
,{ib,numb}];
Zff = Zff<>"\n";
}

a=b1n;pa=1;b=bsr;pb=2;
DC = a^pa*b^pb;
Zff = {};
Zgg = {};

fg = "ZA";

T0 = CmxA*DC;
ff = "ZAmx";
expand;
termsAmx = terms;

T0 = CmyA*DC;
ff = "ZAmy";
expand;
termsAmy = terms;

T0 = CmzA*DC;
ff = "ZAmz";
expand;
termsAmz = terms;

T0 = CzxA*DC;
ff = "ZAzx";
expand;
termsAzx = terms;

T0 = CzyA*DC;
ff = "ZAzy";
expand;
termsAzy = terms;

T0 = CzzA*DC;
ff = "ZAzz";
expand;
termsAzz = terms;

Print[Zgg]
Print[Zff]












