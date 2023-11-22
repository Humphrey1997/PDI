#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
fsLr = S*q*f00/I/T/(ws-ksz*vz-r*wc)*(vz*Esz*BesselJ[r,as*x]+vp/2*((Esx-I*Esy)*Exp[I*des]*BesselJ[r-S,as*x]+(Esx+I*Esy)*Exp[-I*des]*BesselJ[r+S,as*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjr = S*q*f00/I/T/(w0-k0z*vz-(+nr)*wc)*(vz*E0z*BesselJ[+nr,a0*x]+vp/2*((E0x-I*E0y)*Exp[I*de0]*BesselJ[+nr-S,a0*x]+(E0x+I*E0y)*Exp[-I*de0]*BesselJ[+nr+S,a0*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};

f2QLnA0 = 1/I/(w2-k2z*vz-n*wc)*1/2*(dx0*(I*ksx*fsLr)+dy0*(I*ksy*fsLr)+dth0*(I*S*r*fsLr)+dmiu0*D[fsLr,miu]+dpz0*D[fsLr,pz])/.{t->nr};
f2QLnB0 = 1/I/(w2-k2z*vz-n*wc)*1/2*(dxs*(I*k0x*f0Lnjr)+dys*(I*k0y*f0Lnjr)+dths*(I*S*nr*f0Lnjr)+dmius*D[f0Lnjr,miu]+dpzs*D[f0Lnjr,pz])/.{t->r};

f2QLnA1 = f2QLnA0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k2z->n2z,k0z->n0z,ksx->asx*wc/vt,k2x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k2y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k2p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-c2 u+w2-n wc)->b2n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-(+nr) wc)->b0nr}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de2]->Ei2,Exp[-I*de2]->Ei2c};
f2QLnb2 = f2QLnB0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k2z->n2z,k0z->n0z,ksx->asx*wc/vt,k2x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k2y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k2p->a1*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-c2 u+w2-n wc)->b2n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-(+nr) wc)->b0nr}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de2]->Ei2,Exp[-I*de2]->Ei2c};

f2QLnA2 = f2QLnA1/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[+nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4};
f2QLnB2 = f2QLnb2/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[+nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4};

(*special simplify*)
f2QLnA3 = f2QLnA2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};
f2QLnB3 = f2QLnB2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};

JxyA= Simplify[f2QLnA3*m*m*vp/wc*(I/ep0/w2)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzA= Simplify[f2QLnA3*m*m*vp/wc*(I/ep0/w2)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmyA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmzA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]
czxA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
czyA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
czzA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]

JxyB= Simplify[f2QLnB3*m*m*vp/wc*(I/ep0/w2)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzB= Simplify[f2QLnB3*m*m*vp/wc*(I/ep0/w2)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmyB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
cmzB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]
czxB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->1,Esy->0,Esz->0}/.{E0->beta*m*w0*vt/q}]
czyB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->1,Esz->0}/.{E0->beta*m*w0*vt/q}]
czzB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/w2)/.{Esx->0,Esy->0,Esz->1}/.{E0->beta*m*w0*vt/q}]




