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

a=b1n;pa=1;b=bsp;pb=2;c=b1r;pc=3;
Dw = w0*x;
DC = a^pa*b^pb*c^pc*Dw;
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




