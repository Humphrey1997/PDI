#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
f1Lr = S*q*f00/I/T/(w1-k1z*vz-r*wc)*(vz*E1z*BesselJ[r,a1*x]+vp/2*((E1x-I*E1y)*Exp[I*de1]*BesselJ[r-S,a1*x]+(E1x+I*E1y)*Exp[-I*de1]*BesselJ[r+S,a1*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lpjr = S*q*f00/I/T/(w0-k0z*vz-(+pr)*wc)*(vz*E0z*BesselJ[+pr,a0*x]+vp/2*((E0x-I*E0y)*Exp[I*de0]*BesselJ[+pr-S,a0*x]+(E0x+I*E0y)*Exp[-I*de0]*BesselJ[+pr+S,a0*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};

dx1 = dx/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dy1 = dy/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dmiu1 = dmiu/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dth1 = dth/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};
dpz1 = dpz/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,Ejx->E1x,Ejy->E1y,Ejz->E1z,aj->a1,wj->w1,dej->de1};

fsQLpA = 1/I/(ws-ksz*vz-p*wc)*1/2*(dx0*(I*k1x*f1Lr)+dy0*(I*k1y*f1Lr)+dth0*(I*S*r*f1Lr)+dmiu0*D[f1Lr,miu]+dpz0*D[f1Lr,pz])/.{t->pr}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
fsQLpB = 1/I/(ws-ksz*vz-p*wc)*1/2*(dx1*(I*k0x*f0Lpjr)+dy1*(I*k0y*f0Lpjr)+dth1*(I*S*pr*f0Lpjr)+dmiu1*D[f0Lpjr,miu]+dpz1*D[f0Lpjr,pz])/.{t->r}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

f1NLnA0 =  1/I/(w1-k1z*vz-n*wc)*1/2*(dx0c*(I*ksx*fsQLpA)+dy0c*(I*ksy*fsQLpA)+dth0c*(I*S*p*fsQLpA)+dmiu0c*D[fsQLpA,miu]+dpz0c*D[fsQLpA,pz])/.{t->pn};
f1NLnB0 =  1/I/(w1-k1z*vz-n*wc)*1/2*(dx0c*(I*ksx*fsQLpB)+dy0c*(I*ksy*fsQLpB)+dth0c*(I*S*p*fsQLpB)+dmiu0c*D[fsQLpB,miu]+dpz0c*D[fsQLpB,pz])/.{t->pn};

f1NLnA1 = f1NLnA0/.{E^((-(pz^2/(2 m))-miu wc)/T)->1}/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c1 u+w1-n wc)->b1n,(-cs u+ws-p wc)->bsp,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-(+pr) wc)->b0pr}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};
f1NLnB1 = f1NLnB0/.{E^((-(pz^2/(2 m))-miu wc)/T)->1}/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c1 u+w1-n wc)->b1n,(-cs u+ws-p wc)->bsp,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-(+pr) wc)->b0pr}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};

f1NLnA2 = f1NLnA1/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,BesselJ[-3+r,a1 x]->j1r5,BesselJ[3+r,a1 x]->j1r6,
BesselJ[+pr,a0 x]->j0pr0,BesselJ[-1+pr,a0 x]->j0pr1,BesselJ[1+pr,a0 x]->j0pr2,BesselJ[-2+pr,a0 x]->j0pr3,BesselJ[2+pr,a0 x]->j0pr4,BesselJ[-3+pr,a0 x]->j0pr5,BesselJ[3+pr,a0 x]->j0pr6,
BesselJ[+pn,a0 x]->j0pn0,BesselJ[-1+pn,a0 x]->j0pn1,BesselJ[1+pn,a0 x]->j0pn2};

f1NLnB2 = f1NLnB1/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,BesselJ[-3+r,a1 x]->j1r5,BesselJ[3+r,a1 x]->j1r6,
BesselJ[+pr,a0 x]->j0pr0,BesselJ[-1+pr,a0 x]->j0pr1,BesselJ[1+pr,a0 x]->j0pr2,BesselJ[-2+pr,a0 x]->j0pr3,BesselJ[2+pr,a0 x]->j0pr4,BesselJ[-3+pr,a0 x]->j0pr5,BesselJ[3+pr,a0 x]->j0pr6,
BesselJ[+pn,a0 x]->j0pn0,BesselJ[-1+pn,a0 x]->j0pn1,BesselJ[1+pn,a0 x]->j0pn2};

f1NLnA3 = f1NLnA2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};
f1NLnB3 = f1NLnB2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};

JxyA= f1NLnA3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt};
JzzA= f1NLnA3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt};
cmxA = Simplify[JxyA*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmyA = Simplify[JxyA*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmzA = Simplify[JxyA*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]
czxA = Simplify[JzzA*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
czyA = Simplify[JzzA*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
czzA = Simplify[JzzA*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]

JxyB= f1NLnB3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt};
JzzB= f1NLnB3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt};
cmxB = Simplify[JxyB*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmyB = Simplify[JxyB*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
cmzB = Simplify[JxyB*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]
czxB = Simplify[JzzB*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->1,E1y->0,E1z->0}/.{E0->beta*m*w0*vt/q}]
czyB = Simplify[JzzB*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->1,E1z->0}/.{E0->beta*m*w0*vt/q}]
czzB = Simplify[JzzB*Sqrt[Pi]/(beta^2*wp2*w0/w1)/.{E1x->0,E1y->0,E1z->1}/.{E0->beta*m*w0*vt/q}]



