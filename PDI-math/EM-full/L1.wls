#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
f1Ln = S*q*f00/I/T/(w1-k1z*vz-n*wc)*(vz*E1z*BesselJ[n,a1*x]+vp/2*((E1x-I*E1y)*Exp[I*de1]*BesselJ[n-S,a1*x]+(E1x+I*E1y)*Exp[-I*de1]*BesselJ[n+S,a1*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

f1Ln1 = f1Ln/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};

f1Ln2 = f1Ln1/.{BesselJ[n,a1 x]->j1n0,BesselJ[-1+n,a1 x]->j1n1,BesselJ[1+n,a1 x]->j1n2}/.
{(-n1z u vt+w1-n wc)->b1n};

JxyA= Simplify[f1Ln2*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzA= Simplify[f1Ln2*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(wp2/w1)/.{E1x->1,E1y->0,E1z->0}]
cmyA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(wp2/w1)/.{E1x->0,E1y->1,E1z->0}]
cmzA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(wp2/w1)/.{E1x->0,E1y->0,E1z->1}]
czxA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(wp2/w1)/.{E1x->1,E1y->0,E1z->0}]
czyA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(wp2/w1)/.{E1x->0,E1y->1,E1z->0}]
czzA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(wp2/w1)/.{E1x->0,E1y->0,E1z->1}]



