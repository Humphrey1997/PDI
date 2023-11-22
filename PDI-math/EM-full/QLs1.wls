#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
f1Lr = S*q*f00/I/T/(w1-k1z*vz-r*wc)*(vz*E1z*BesselJ[r,a1*x]+vp/2*((E1x-I*E1y)*Exp[I*de1]*BesselJ[r-S,a1*x]+(E1x+I*E1y)*Exp[-I*de1]*BesselJ[r+S,a1*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjr = S*q*f00/I/T/(w0-k0z*vz-(+nr)*wc)*(vz*E0z*BesselJ[+nr,a0*x]+vp/2*((E0x-I*E0y)*Exp[I*de0]*BesselJ[+nr-S,a0*x]+(E0x+I*E0y)*Exp[-I*de0]*BesselJ[+nr+S,a0*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};

dx1 = dx/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dy1 = dy/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dmiu1 = dmiu/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dth1 = dth/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dpz1 = dpz/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};

fs1QLnA0 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1Lr)+dy0*(I*k1y*f1Lr)+dth0*(I*S*r*f1Lr)+dmiu0*D[f1Lr,miu]+dpz0*D[f1Lr,pz])/.{t->+nr};
fs1QLnB0 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx1*(I*k0x*f0Lnjr)+dy1*(I*k0y*f0Lnjr)+dth1*(I*S*nr*f0Lnjr)+dmiu1*D[f0Lnjr,miu]+dpz1*D[f0Lnjr,pz])/.{t->r};

fs1QLnA1 = fs1QLnA0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-(+nr) wc)->b0nr}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};
fs1QLnB1 = fs1QLnB0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-(+nr) wc)->b0nr}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};

fs1QLnA2 = fs1QLnA1/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,
BesselJ[+nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4};
fs1QLnB2 = fs1QLnB1/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,
BesselJ[+nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4};

(*special simplify*)
fs1QLnA3 = fs1QLnA2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};
fs1QLnB3 = fs1QLnB2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};

JxyA= Simplify[fs1QLnA3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzA= Simplify[fs1QLnA3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmyA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmzA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]
czxA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
czyA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
czzA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]

JxyB= Simplify[fs1QLnB3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzB= Simplify[fs1QLnB3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmyB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmzB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]
czxB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
czyB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
czzB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]

