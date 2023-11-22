#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};

fsLr = 1/I/(ws-ksz*vz-t*wc)*(dmius*D[f00,miu]+dpzs*D[f00,pz])/.{t->r}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjrc = -1/I/(w0-k0z*vz-t*wc)*(dmiu0c*D[f00,miu]+dpz0c*D[f00,pz])/.{t->rn}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

f1QLn1 = 1/I/(w1-k1z*vz-n*wc)*1/2*(dx0c*(I*ksx*fsLr)+dy0c*(I*ksy*fsLr)+dth0c*(I*S*r*fsLr)+dmiu0c*D[fsLr,miu]+dpz0c*D[fsLr,pz])/.{t->rn};
f1QLn2 = 1/I/(w1-k1z*vz-n*wc)*1/2*(dxs*(-I*k0x*f0Lnjrc)+dys*(-I*k0y*f0Lnjrc)+dths*(-I*S*rn*f0Lnjrc)+dmius*D[f0Lnjrc,miu]+dpzs*D[f0Lnjrc,pz])/.{t->r};

f1QLn1A = f1QLn1/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
f1QLn2A = f1QLn2/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

f1QLn1B=f1QLn1A/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c1 u+w1-n wc)->b1n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-rn wc)->b0rn};
f1QLn2B=f1QLn2A/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c1 u+w1-n wc)->b1n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-rn wc)->b0rn};

(*special simplify*)
f1QLn1C = f1QLn1B/.{a0x->a0,a0y->0};
f1QLn2C = f1QLn2B/.{a0x->a0,a0y->0};

NA1 = f1QLn1C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a1*x]*(S*q)/ep0/k1^2/.{ep0->n0*q*q/m/wp2}/.{k1->n1}/.{BesselJ[n,a1*x]->j1n0};
NA2 = f1QLn2C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a1*x]*(S*q)/ep0/k1^2/.{ep0->n0*q*q/m/wp2}/.{k1->n1}/.{BesselJ[n,a1*x]->j1n0};
NB1 = NA1*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2 = NA2*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1 = Simplify[NB1/(wp2*yita0/n1^2/vt^2)]
NC2 = Simplify[NB2/(wp2*yita0/n1^2/vt^2)]














