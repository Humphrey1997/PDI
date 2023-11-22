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

a=b2n;pa=1;b=b0nr;pb=2;
DC = a^pa*b^pb;
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






