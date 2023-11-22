#!/usr/bin/env wolframscript
(* ::Package:: *)

a=b2n;pa=1;b=bsp;pb=2;c=b2r;pc=3;
DC1 = a^pa*b^pb*c^pc;
N0 = NC1*DC1;
count = 0;
termsA = {};
ZfA = {};
ca = CoefficientList[N0,a];
numa = Part[Dimensions[ca],1];
Do[
TA = Part[ca,ia];
cb = CoefficientList[TA,b];
numb = Part[Dimensions[cb],1];
Do[
TB = Part[cb,ib];
cc = CoefficientList[TB,c];
numc = Part[Dimensions[cc],1];
Do[
TC = Part[cc,ic];
cu = FullSimplify[CoefficientList[TC,u]];
numu = Part[Dimensions[cu],1];
Do[
term = Part[cu,iu];
If[ToString[term]!="0",
indexa = StringPart[ToString[a],2;;-1];
indexpa = ToString[pa-ia+1];
indexb = StringPart[ToString[b],2;;-1];
indexpb = ToString[pb-ib+1];
indexc = StringPart[ToString[c],2;;-1];
indexpc = ToString[pc-ic+1];
indexu = ToString[iu-1];
Zf = "ZA"<>ToString[count+1];
count = count + 1;
termsA = Append[termsA,term];
G = "g3("<>indexu<>",zeta"<>indexa<>","<>indexpa<>",zeta"<>indexb<>","<>indexpb<>",zeta"<>indexc<>","<>indexpc<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa<>"/c"<>StringPart[ToString[b],2]<>"^"<>indexpb<>"/c"<>StringPart[ToString[c],2]<>"^"<>indexpc;
ZfA = ZfA<>Zf<>"="<>G<>CS<>";\n";
,Print["Jump"]];
,{iu,numu}];
,{ic,numc}];
,{ib,numb}];
,{ia,numa}];
Print[ZfA];

a=b2n;pa=1;b=bsp;pb=2;c=b0rp;pc=3;
DC2 = a^pa*b^pb*c^pc;
N0 = NC2*DC2;
count = 0;
termsB = {};
ZfB = {};
ca = CoefficientList[N0,a];
numa = Part[Dimensions[ca],1];
Do[
TA = Part[ca,ia];
cb = CoefficientList[TA,b];
numb = Part[Dimensions[cb],1];
Do[
TB = Part[cb,ib];
cc = CoefficientList[TB,c];
numc = Part[Dimensions[cc],1];
Do[
TC = Part[cc,ic];
cu = FullSimplify[CoefficientList[TC,u]];
numu = Part[Dimensions[cu],1];
Do[
term = Part[cu,iu];
If[ToString[term]!="0",
indexa = StringPart[ToString[a],2;;-1];
indexpa = ToString[pa-ia+1];
indexb = StringPart[ToString[b],2;;-1];
indexpb = ToString[pb-ib+1];
indexc = StringPart[ToString[c],2;;-1];
indexpc = ToString[pc-ic+1];
indexu = ToString[iu-1];
Zf = "ZB"<>ToString[count+1];
count = count + 1;
termsB = Append[termsB,term];
G = "g3("<>indexu<>",zeta"<>indexa<>","<>indexpa<>",zeta"<>indexb<>","<>indexpb<>",zeta"<>indexc<>","<>indexpc<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa<>"/c"<>StringPart[ToString[b],2]<>"^"<>indexpb<>"/c"<>StringPart[ToString[c],2]<>"^"<>indexpc;
ZfB = ZfB<>Zf<>"="<>G<>CS<>";\n";
,Print["Jump"]];
,{iu,numu}];
,{ic,numc}];
,{ib,numb}];
,{ia,numa}];
Print[ZfB];




