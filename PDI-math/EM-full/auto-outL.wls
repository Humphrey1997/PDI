#!/usr/bin/env wolframscript
(* ::Package:: *)

TMA :={
terms = ToExpression["termsA"<>dir];
yA = "Y"<>"A"<>dir;
zA = "Z"<>"A"<>dir;
YfA = {};
FfA="";
numt = Part[Dimensions[terms],1];
Do[
Mterm = ToMatlab[Part[terms,it],1000000];
Mterm = StringReplace[Mterm,"((-1).*"->"(-"];
Mterm = StringReplace[Mterm,"+(-1).*"->"-"];
Mterm = StringReplace[Mterm,"1i*(-1)"->"-1i"];
Mterm = StringReplace[Mterm,"1i*(-2)"->"-2i"];
Mterm = StringReplace[Mterm,"1i*2"->"2i"];
Mterm = StringReplace[Mterm,"1i*(-4)"->"-4i"];
Mterm = StringReplace[Mterm,"1i*4"->"4i"];
Mterm = StringReplace[Mterm,"(-1)+"->"-1+"];
Mterm = StringReplace[Mterm,"(-1)-"->"-1-"];
Mterm = StringReplace[Mterm,"(-2)+"->"-2+"];
Mterm = StringReplace[Mterm,"(-2)-"->"-2-"];
YfA = YfA<>yA<>ToString[it]<>" = "<>Mterm;
FfA = FfA<>yA<>ToString[it]<>"*"<>zA<>ToString[it]<>"+";
,{it,numt}]
}


TM:={
TMA;
WriteString[fo,YfA];
Ff = StringDrop[FfA,-1];
Fall = Fall<>"c"<>dir<>" = exp(-x.^2).*("<>Ff<>");\n";
}

fo=OpenWrite["D:\\Documents\\Wolfram Mathematica\\2022\\1PDI -full\\outL.txt"];
Fall = "";
dir = "mx";
TM;
dir = "my";
TM;
dir = "mz";
TM;
dir = "zx";
TM;
dir = "zy";
TM;
dir = "zz";
TM;

WriteString[fo,"\n"];
WriteString[fo,Fall];
Close[fo]


