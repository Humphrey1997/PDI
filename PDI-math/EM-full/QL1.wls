#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
fsLr = S*q*f00/I/T/(ws-ksz*vz-r*wc)*(vz*Esz*BesselJ[r,as*x]+vp/2*((Esx-I*Esy)*Exp[I*des]*BesselJ[r-S,as*x]+(Esx+I*Esy)*Exp[-I*des]*BesselJ[r+S,as*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjrc = S*q*f00/(-I)/T/(w0-k0z*vz-(+rn)*wc)*(vz*E0zc*BesselJ[+rn,a0*x]+vp/2*((E0xc+I*E0yc)*Exp[-I*de0]*BesselJ[+rn-S,a0*x]+(E0xc-I*E0yc)*Exp[I*de0]*BesselJ[+rn+S,a0*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};

f1QLnA0 = 1/I/(w1-k1z*vz-n*wc)*1/2*(dx0c*(I*ksx*fsLr)+dy0c*(I*ksy*fsLr)+dth0c*(I*S*r*fsLr)+dmiu0c*D[fsLr,miu]+dpz0c*D[fsLr,pz])/.{t->+rn};
f1QLnB0 = 1/I/(w1-k1z*vz-n*wc)*1/2*(dxs*(-I*k0x*f0Lnjrc)+dys*(-I*k0y*f0Lnjrc)+dths*(-I*S*rn*f0Lnjrc)+dmius*D[f0Lnjrc,miu]+dpzs*D[f0Lnjrc,pz])/.{t->r};

f1QLnA1 = f1QLnA0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c1 u+w1-n wc)->b1n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-(+rn) wc)->b0rn}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};
f1QLnB1 = f1QLnB0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k0z->n0z,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c1 u+w1-n wc)->b1n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-(+rn) wc)->b0rn}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c};

f1QLnA2 = f1QLnA1/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[+rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4};
f1QLnB2 = f1QLnB1/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[+rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4};

(*special simplify*)
f1QLnA3 = f1QLnA2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};
f1QLnB3 = f1QLnB2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};

JxyA= Simplify[f1QLnA3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzA= Simplify[f1QLnA3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmyA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmzA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]
czxA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
czyA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
czzA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]

JxyB= Simplify[f1QLnB3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzB= Simplify[f1QLnB3*m*m*vp/wc*(I/ep0/w1)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmyB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmzB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]
czxB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
czyB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
czzB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w1)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]




