#!/usr/bin/env wolframscript
(* ::Package:: *)

terms = termsA;
YfA = {};
numt = Part[Dimensions[terms],1];
Do[
Mterm = ToMatlab[Part[terms,it],1000000];
Mterm = StringReplace[Mterm,"((-1).*"->"(-"];
Mterm = StringReplace[Mterm,"+(-1).*"->"-"];
Mterm = StringReplace[Mterm,"1i*(-1)"->"-1i"];
Mterm = StringReplace[Mterm,"1i*(-2)"->"-2i"];
Mterm = StringReplace[Mterm,"1i*(-4)"->"-4i"];
Mterm = StringReplace[Mterm,"(-1)+"->"-1+"];
Mterm = StringReplace[Mterm,"(-1)-"->"-1-"];
Mterm = StringReplace[Mterm,"(-2)+"->"-2+"];
Mterm = StringReplace[Mterm,"(-2)-"->"-2-"];
YfA = YfA<>"YA"<>ToString[it]<>" = "<>Mterm;
,{it,numt}]


FfA="";
Do[
FfA = FfA<>"YA"<>ToString[it]<>"*"<>"ZA"<>ToString[it]<>"+";
,{it,numt}]
Print[FfA];
fo=OpenWrite["D:\\Documents\\Wolfram Mathematica\\2022\\2PDI - full\\outL.txt"]
WriteString[fo,YfA]
Close[fo]

















