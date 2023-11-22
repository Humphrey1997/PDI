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
Mterm = StringReplace[Mterm,"(1i*2)"->"2i"];
Mterm = StringReplace[Mterm,"1i*(-4)"->"-4i"];
Mterm = StringReplace[Mterm,"(1i*4)"->"4i"];
Mterm = StringReplace[Mterm,"(-1)+"->"-1+"];
Mterm = StringReplace[Mterm,"(-1)-"->"-1-"];
Mterm = StringReplace[Mterm,"(-2)+"->"-2+"];
Mterm = StringReplace[Mterm,"(-2)-"->"-2-"];
YfA = YfA<>"YA"<>ToString[it]<>" = "<>Mterm;
,{it,numt}]

terms = termsB;
YfB = {};
numt = Part[Dimensions[terms],1];
Do[
Mterm = ToMatlab[Part[terms,it],1000000];
Mterm = StringReplace[Mterm,"((-1).*"->"(-"];
Mterm = StringReplace[Mterm,"+(-1).*"->"-"];
Mterm = StringReplace[Mterm,"1i*(-1)"->"-1i"];
Mterm = StringReplace[Mterm,"1i*(-2)"->"-2i"];
Mterm = StringReplace[Mterm,"(1i*2)"->"2i"];
Mterm = StringReplace[Mterm,"1i*(-4)"->"-4i"];
Mterm = StringReplace[Mterm,"(1i*4)"->"4i"];
Mterm = StringReplace[Mterm,"(-1)+"->"-1+"];
Mterm = StringReplace[Mterm,"(-1)-"->"-1-"];
Mterm = StringReplace[Mterm,"(-2)+"->"-2+"];
Mterm = StringReplace[Mterm,"(-2)-"->"-2-"];
YfB = YfB<>"YB"<>ToString[it]<>" = "<>Mterm;
,{it,numt}]

FfA="";
FfB="";
Do[
FfA = FfA<>"YA"<>ToString[it]<>"*"<>"ZA"<>ToString[it]<>"+";
FfB = FfB<>"YB"<>ToString[it]<>"*"<>"ZB"<>ToString[it]<>"+";
,{it,numt}]
Print[FfA<>FfB];
fo=OpenWrite["D:\\Code\\Mathematica\\202209\\ES-full\\out.txt"]
WriteString[fo,YfA]
WriteString[fo,YfB]
Close[fo]









