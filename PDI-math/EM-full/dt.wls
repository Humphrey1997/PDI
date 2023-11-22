#!/usr/bin/env wolframscript
(* ::Package:: *)

(*for electron S=1, for ion S=-1*)
S = 1;

dx = q/m/wc/wj*((wj-kjz*vz)*Ejy*BesselJ[t,aj*x]+kjy*vz*Ejz*BesselJ[t,aj*x]+vp/2*(kjy*Ejx-kjx*Ejy)*(Exp[I*dej]*BesselJ[t-S,aj*x]+Exp[-I*dej]*BesselJ[t+S,aj*x]));
dy = -q/m/wc/wj*((wj-kjz*vz)*Ejx*BesselJ[t,aj*x]+kjx*vz*Ejz*BesselJ[t,aj*x]+I*vp/2*(kjy*Ejx-kjx*Ejy)*(Exp[I*dej]*BesselJ[t-S,aj*x]-Exp[-I*dej]*BesselJ[t+S,aj*x]));
dmiu = -S*q/wc/wj*(t*wc*vz*Ejz*BesselJ[t,aj*x]+vp/2*(wj-kjz*vz)*((Ejx-I*Ejy)*Exp[I*dej]*BesselJ[t-S,aj*x]+(Ejx+I*Ejy)*Exp[-I*dej]*BesselJ[t+S,aj*x]));
dth = -S*q/m/vp/wj*(I/2*(wj-kjz*vz)*((Ejx-I*Ejy)*Exp[I*dej]*BesselJ[t-S,aj*x]-(Ejx+I*Ejy)*Exp[-I*dej]*BesselJ[t+S,aj*x])+I*S*kjp*vz*Ejz*(BesselJ[t-1,aj*x]-BesselJ[t+1,aj*x])/2+vp*(kjy*Ejx-kjx*Ejy)*BesselJ[t,aj*x]);
dpz = -S*q/wj*((wj-t*wc)*Ejz*BesselJ[t,aj*x]+kjz*vp/2*((Ejx-I*Ejy)*Exp[I*dej]*BesselJ[t-S,aj*x]+(Ejx+I*Ejy)*Exp[-I*dej]*BesselJ[t+S,aj*x]));


dxc = q/m/wc/wj*((wj-kjz*vz)*Ejyc*BesselJ[t,aj*x]+kjy*vz*Ejzc*BesselJ[t,aj*x]+vp/2*(kjy*Ejxc-kjx*Ejyc)*(Exp[-I*dej]*BesselJ[t-S,aj*x]+Exp[I*dej]*BesselJ[t+S,aj*x]));
dyc = -q/m/wc/wj*((wj-kjz*vz)*Ejxc*BesselJ[t,aj*x]+kjx*vz*Ejzc*BesselJ[t,aj*x]-I*vp/2*(kjy*Ejxc-kjx*Ejyc)*(Exp[-I*dej]*BesselJ[t-S,aj*x]-Exp[I*dej]*BesselJ[t+S,aj*x]));
dmiuc = -S*q/wc/wj*(t*wc*vz*Ejzc*BesselJ[t,aj*x]+vp/2*(wj-kjz*vz)*((Ejxc+I*Ejyc)*Exp[-I*dej]*BesselJ[t-S,aj*x]+(Ejxc-I*Ejyc)*Exp[I*dej]*BesselJ[t+S,aj*x]));
dthc = -S*q/m/vp/wj*(-I/2*(wj-kjz*vz)*((Ejxc+I*Ejyc)*Exp[-I*dej]*BesselJ[t-S,aj*x]-(Ejxc-I*Ejyc)*Exp[I*dej]*BesselJ[t+S,aj*x])-I*S*kjp*vz*Ejzc*(BesselJ[t-1,aj*x]-BesselJ[t+1,aj*x])/2+vp*(kjy*Ejxc-kjx*Ejyc)*BesselJ[t,aj*x]);
dpzc = -S*q/wj*((wj-t*wc)*Ejzc*BesselJ[t,aj*x]+kjz*vp/2*((Ejxc+I*Ejyc)*Exp[-I*dej]*BesselJ[t-S,aj*x]+(Ejxc-I*Ejyc)*Exp[I*dej]*BesselJ[t+S,aj*x]));
