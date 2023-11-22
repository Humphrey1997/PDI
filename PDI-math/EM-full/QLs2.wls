#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
f2Lr = S*q*f00/I/T/(w2-k2z*vz-r*wc)*(vz*E2z*BesselJ[r,a2*x]+vp/2*((E2x-I*E2y)*Exp[I*de2]*BesselJ[r-S,a2*x]+(E2x+I*E2y)*Exp[-I*de2]*BesselJ[r+S,a2*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjrc = S*q*f00/(-I)/T/(w0-k0z*vz-(+rn)*wc)*(vz*E0zc*BesselJ[+rn,a0*x]+vp/2*((E0xc+I*E0yc)*Exp[-I*de0]*BesselJ[+rn-S,a0*x]+(E0xc-I*E0yc)*Exp[I*de0]*BesselJ[+rn+S,a0*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};

dx2 = dx/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,Ejx->E2x,Ejy->E2y,Ejz->E2z,aj->a2,wj->w2,dej->de2};
dy2 = dy/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,Ejx->E2x,Ejy->E2y,Ejz->E2z,aj->a2,wj->w2,dej->de2};
dmiu2 = dmiu/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,Ejx->E2x,Ejy->E2y,Ejz->E2z,aj->a2,wj->w2,dej->de2};
dth2 = dth/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,Ejx->E2x,Ejy->E2y,Ejz->E2z,aj->a2,wj->w2,dej->de2};
dpz2 = dpz/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,Ejx->E2x,Ejy->E2y,Ejz->E2z,aj->a2,wj->w2,dej->de2};

fs2QLnA0 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0c*(I*k2x*f2Lr)+dy0c*(I*k2y*f2Lr)+dth0c*(I*S*r*f2Lr)+dmiu0c*D[f2Lr,miu]+dpz0c*D[f2Lr,pz])/.{t->+rn};
fs2QLnB0 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx2*(-I*k0x*f0Lnjrc)+dy2*(-I*k0y*f0Lnjrc)+dth2*(-I*S*rn*f0Lnjrc)+dmiu2*D[f0Lnjrc,miu]+dpz2*D[f0Lnjrc,pz])/.{t->r};

fs2QLnA1 = fs2QLnA0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k2z->n2z,k0z->n0z,ksx->asx*wc/vt,k2x->a2x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k2y->a2y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k2p->a2*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c2 u+w2-r wc)->b2r,(-c0 u+w0-(+rn) wc)->b0rn}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de2]->Ei2,Exp[-I*de2]->Ei2c};
fs2QLnB1 = fs2QLnB0/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k2z->n2z,k0z->n0z,ksx->asx*wc/vt,k2x->a2x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k2y->a2y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k2p->a2*wc/vt,k0p->a0*wc/vt}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c2 u+w2-r wc)->b2r,(-c0 u+w0-(+rn) wc)->b0rn}/.
{E0x+I*E0y->e0p*E0,E0x-I*E0y->e0s*E0,E0xc+I*E0yc->e0sc*E0,E0xc-I*E0yc->e0pc*E0}/.
{E0z->e0z*E0,E0y->e0y*E0,E0x->e0x*E0,E0zc->e0zc*E0,E0yc->e0yc*E0,E0xc->e0xc*E0}/.
{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de2]->Ei2,Exp[-I*de2]->Ei2c};

fs2QLnA2 = fs2QLnA1/.{BesselJ[r,a2 x]->j2r0,BesselJ[-1+r,a2 x]->j2r1,BesselJ[1+r,a2 x]->j2r2,BesselJ[-2+r,a2 x]->j2r3,BesselJ[2+r,a2 x]->j2r4,
BesselJ[+rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4};
fs2QLnB2 = fs2QLnB1/.{BesselJ[r,a2 x]->j2r0,BesselJ[-1+r,a2 x]->j2r1,BesselJ[1+r,a2 x]->j2r2,BesselJ[-2+r,a2 x]->j2r3,BesselJ[2+r,a2 x]->j2r4,
BesselJ[+rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4};

(*special simplify*)
fs2QLnA3 = fs2QLnA2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};
fs2QLnB3 = fs2QLnB2/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1};

JxyA= Simplify[fs2QLnA3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzA= Simplify[fs2QLnA3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->1,E2y->0,E2z->0}/.{E0->beta*m*w0*vt/q}]
cmyA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->1,E2z->0}/.{E0->beta*m*w0*vt/q}]
cmzA = Simplify[JxyA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->0,E2z->1}/.{E0->beta*m*w0*vt/q}]
czxA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->1,E2y->0,E2z->0}/.{E0->beta*m*w0*vt/q}]
czyA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->1,E2z->0}/.{E0->beta*m*w0*vt/q}]
czzA = Simplify[JzzA*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->0,E2z->1}/.{E0->beta*m*w0*vt/q}]

JxyB= Simplify[fs2QLnB3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vp*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
JzzB= Simplify[fs2QLnB3*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt,vz->u*vt}];
cmxB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->1,E2y->0,E2z->0}/.{E0->beta*m*w0*vt/q}]
cmyB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->1,E2z->0}/.{E0->beta*m*w0*vt/q}]
cmzB = Simplify[JxyB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->0,E2z->1}/.{E0->beta*m*w0*vt/q}]
czxB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->1,E2y->0,E2z->0}/.{E0->beta*m*w0*vt/q}]
czyB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->1,E2z->0}/.{E0->beta*m*w0*vt/q}]
czzB = Simplify[JzzB*Sqrt[Pi]*Exp[u^2+x^2]/(beta*wp2*w0/ws)/.{E2x->0,E2y->0,E2z->1}/.{E0->beta*m*w0*vt/q}]




