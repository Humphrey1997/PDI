#!/usr/bin/env wolframscript
(* ::Package:: *)

expand :={
terms = {};
count = 0;
Ta = T0;
indexa = StringPart[ToString[a],2;;-1];
indexb = StringPart[ToString[b],2;;-1];
indexc = StringPart[ToString[c],2;;-1];
indexpa = "1";
cb = CoefficientList[Ta,b];
numb = Part[Dimensions[cb],1];
Do[
Tb = Part[cb,ib];
cc = CoefficientList[Tb,c];
numc = Part[Dimensions[cc],1];
Do[
Tc = Part[cc,ic];
cu = CoefficientList[Tc,u];
numu = Part[Dimensions[cu],1];
Do[
term = Simplify[Part[cu,iu]]/Dw;
If[ToString[term]!="0",
indexpb = ToString[pb-ib+1];
indexpc = ToString[pc-ic+1];
indexpu = ToString[iu-1];
Zg = fg<>indexpu<>indexpa<>indexpb<>indexpc;
Zf = ff<>ToString[count+1];
count = count + 1;
G = "g3("<>indexpu<>",zeta"<>indexa<>","<>indexpa<>",zeta"<>indexb<>","<>indexpb<>",zeta"<>indexc<>","<>indexpc<>")";
CS = "/c"<>StringPart[ToString[a],2]<>"^"<>indexpa<>"/c"<>StringPart[ToString[b],2]<>"^"<>indexpb<>"/c"<>StringPart[ToString[c],2]<>"^"<>indexpc;
Zff = Zff<>Zf<>"="<>Zg<>";";
Zgg = Zgg<>Zg<>"="<>G<>CS<>";\n";
terms = Append[terms,term];
Print[Zg];
,Print["Jump"]];
,{iu,numu}];
,{ic,numc}];
,{ib,numb}];
Zff = Zff<>"\n";
Print[ff];
Print[DateString[]];
}

a=b1n;pa=1;b=bsp;pb=2;c=b0pr;pc=3;
Dw = w1*x;
DC = a^pa*b^pb*c^pc*Dw;
Zff = {};
Zgg = {};

fg = "ZB";

T0 = CmxB*DC;
ff = "ZBmx";
expand;
termsBmx = terms;

T0 = CmyB*DC;
ff = "ZBmy";
expand;
termsBmy = terms;

T0 = CmzB*DC;
ff = "ZBmz";
expand;
termsBmz = terms;

T0 = CzxB*DC;
ff = "ZBzx";
expand;
termsBzx = terms;

T0 = CzyB*DC;
ff = "ZBzy";
expand;
termsBzy = terms;

T0 = CzzB*DC;
ff = "ZBzz";
expand;
termsBzz = terms;

Print[Zgg]
Print[Zff]

