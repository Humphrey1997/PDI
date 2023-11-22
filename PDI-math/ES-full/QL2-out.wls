#!/usr/bin/env wolframscript
(* ::Package:: *)

a=b2n;pa=1;b=bsr;pb=2;
DC1 = a^pa*b^pb;
N0 = NC1*DC1;
Zff = "";
count = 0;
termsA = {};
ca = CoefficientList[N0,a];
numa = Part[Dimensions[ca],1];
Do[
TA = Part[ca,ia];
cb = CoefficientList[TA,b];
numb = Part[Dimensions[cb],1];
Do[
TB = Part[cb,ib];
cu = Simplify[CoefficientList[TB,u]];
numu = Part[Dimensions[cu],1];
Do[
term = Part[cu,iu];
If[ToString[term]!="0",
indexa = StringPart[ToString[a],2;;-1];
indexpa = ToString[pa-ia+1];
indexb = StringPart[ToString[b],2;;-1];
indexpb = ToString[pb-ib+1];
indexu = ToString[iu-1];
Zf = "ZA"<>ToString[count+1];
count = count + 1;
termsA = Append[termsA,term];
G = "g2("<>indexu<>",zeta"<>indexa<>","<>indexpa<>",zeta"<>indexb<>","<>indexpb<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa<>"/c"<>StringPart[ToString[b],2]<>"^"<>indexpb;
Zff = Zff<>Zf<>"="<>G<>CS<>";\n";
,Print["Jump"]];
,{iu,numu}];
,{ib,numb}];
,{ia,numa}];
Print[Zff]

terms2 = 0;
a=b2n;pa=1;b=b0nr;pb=2;
DC2 = a^pa*b^pb;
N0 = NC2*DC2;
Zff = "";
count = 0;
termsB = {};
ca = CoefficientList[N0,a];
numa = Part[Dimensions[ca],1];
Do[
TA = Part[ca,ia];
cb = CoefficientList[TA,b];
numb = Part[Dimensions[cb],1];
Do[
TB = Part[cb,ib];
cu = Simplify[CoefficientList[TB,u]];
numu = Part[Dimensions[cu],1];
Do[
term = Part[cu,iu];
If[ToString[term]!="0",
indexa = StringPart[ToString[a],2;;-1];
indexpa = ToString[pa-ia+1];
indexb = StringPart[ToString[b],2;;-1];
indexpb = ToString[pb-ib+1];
indexu = ToString[iu-1];
Zf = "ZB"<>ToString[count+1];
count = count + 1;
termsB = Append[termsB,term];
G = "g2("<>indexu<>",zeta"<>indexa<>","<>indexpa<>",zeta"<>indexb<>","<>indexpb<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa<>"/c"<>StringPart[ToString[b],2]<>"^"<>indexpb;
Zff = Zff<>Zf<>"="<>G<>CS<>";\n";
,Print["Jump"]];
,{iu,numu}];
,{ib,numb}];
,{ia,numa}];
Print[Zff]




