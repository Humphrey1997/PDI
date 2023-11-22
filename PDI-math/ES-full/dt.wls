#!/usr/bin/env wolframscript
(* ::Package:: *)

(*for electron S=1, for ion S=-1*)
S = 1;

dx = -I*q/m/wc*kjy*BesselJ[t,aj*x]*phj;
dxc = I*q/m/wc*kjy*BesselJ[t,aj*x]*phj;

dy = I*q/m/wc*kjx*BesselJ[t,aj*x]*phj;
dyc = -I*q/m/wc*kjx*BesselJ[t,aj*x]*phj;

dmiu = I*S*q*t*BesselJ[t,aj*x]*phj;
dmiuc = -I*S*q*t*BesselJ[t,aj*x]*phj;

dth = -S*q/2/m/vp*kjp*(BesselJ[t-S,aj*x]-BesselJ[t+S,aj*x])*phj;
dthc = -S*q/2/m/vp*kjp*(BesselJ[t-S,aj*x]-BesselJ[t+S,aj*x])*phj;

dpz = I*S*q*kjz*BesselJ[t,aj*x]*phj;
dpzc = -I*S*q*kjz*BesselJ[t,aj*x]*phj;
