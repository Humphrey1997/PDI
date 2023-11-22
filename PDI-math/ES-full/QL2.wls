#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};

fsLr = 1/I/(ws-ksz*vz-t*wc)*(dmius*D[f00,miu]+dpzs*D[f00,pz])/.{t->r}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjr = 1/I/(w0-k0z*vz-t*wc)*(dmiu0*D[f00,miu]+dpz0*D[f00,pz])/.{t->nr}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

f2QLn1 = 1/I/(w2-k2z*vz-n*wc)*1/2*(dx0*(I*ksx*fsLr)+dy0*(I*ksy*fsLr)+dth0*(I*S*r*fsLr)+dmiu0*D[fsLr,miu]+dpz0*D[fsLr,pz])/.{t->nr};
f2QLn2 = 1/I/(w2-k2z*vz-n*wc)*1/2*(dxs*(I*k0x*f0Lnjr)+dys*(I*k0y*f0Lnjr)+dths*(I*S*nr*f0Lnjr)+dmius*D[f0Lnjr,miu]+dpzs*D[f0Lnjr,pz])/.{t->r};

f2QLn1A = f2QLn1/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
f2QLn2A = f2QLn2/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

f2QLn1B=f2QLn1A/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-c2 u+w2-n wc)->b2n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-(nr) wc)->b0nr};
f2QLn2B=f2QLn2A/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,
BesselJ[nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-c2 u+w2-n wc)->b2n,(-cs u+ws-r wc)->bsr,(-c0 u+w0-(nr) wc)->b0nr};

(*special simplify*)
f2QLn1C = f2QLn1B/.{a0x->a0,a0y->0};
f2QLn2C = f2QLn2B/.{a0x->a0,a0y->0};

NA1 = f2QLn1C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a2*x]*(S*q)/ep0/k2^2/.{ep0->n0*q*q/m/wp2}/.{k2->n2}/.{BesselJ[n,a2*x]->j2n0};
NA2 = f2QLn2C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a2*x]*(S*q)/ep0/k2^2/.{ep0->n0*q*q/m/wp2}/.{k2->n2}/.{BesselJ[n,a2*x]->j2n0};
NB1=NA1*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2=NA2*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1 = Simplify[NB1/(wp2*yita0/n2^2/vt^2)]
NC2 = Simplify[NB2/(wp2*yita0/n2^2/vt^2)]














