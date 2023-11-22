#!/usr/bin/env wolframscript
(* ::Package:: *)

S0 = Import["D:\\Documents\\Wolfram Mathematica\\2022\\1PDI -full\\outT.txt"];
TAB = Table[False,{ix,numx},{ia,numa},{ib,numb},{ic,numc}];
Do[
AS = ToString[ia-1];
Do[
BS = ToString[ib-1];
Do[
CS = ToString[ic-1];
Do[
XS = ToString[ix-1];
SS = "I"<>XS<>AS<>BS<>CS;
TAB[[ix,ia,ib,ic]] = StringContainsQ[S0,SS];
,{ix,numx}]
,{ic,numc}]
,{ib,numb}]
,{ia,numa}]

stA = StringPart[AT,3];
stB = StringPart[BT,3];
stC1 = StringPart[CT,3];
stC2 = StringPart[CT,4];
stC = stC1<>"-"<>stC2;
skA = "a"<>StringPart[AT,2];
skB = "a"<>StringPart[BT,2];
skC = "a"<>StringPart[CT,2];
fi=OpenWrite["D:\\Documents\\Wolfram Mathematica\\2022\\1PDI -full\\outi.txt"];
Do[
XS = ToString[ix-1];
Do[
AS = ToString[ia-1];
Do[
BS = ToString[ib-1];
Do[
CS = ToString[ic-1];
SS = "I"<>XS<>AS<>BS<>CS;
If[TAB[[ix,ia,ib,ic]],
If[Mod[ia-1,2]==0,pA=stA<>"+"<>ToString[(ia-1)/2],pA=stA<>"-"<>ToString[ia/2]];
If[Mod[ib-1,2]==0,pB=stB<>"+"<>ToString[(ib-1)/2],pB=stB<>"-"<>ToString[ib/2]];
If[Mod[ic-1,2]==0,pC=stC<>"+"<>ToString[(ic-1)/2],pC=stC<>"-"<>ToString[ic/2]];
IS = SS<>" = I3("<>XS<>","<>skA<>","<>pA<>","<>skB<>","<>pB<>","<>skC<>","<>pC<>");";
WriteString[fi,IS];
,kkk = 0;];
,{ic,numc}]
,{ib,numb}]
,{ia,numa}]
WriteString[fi,"\n"];
,{ix,numx}]
Close[fi]
