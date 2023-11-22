#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

dx1 = dx/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dy1 = dy/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dmiu1 = dmiu/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dth1 = dth/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dpz1 = dpz/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};

f1Lr = 1/I/(w1-k1z*vz-t*wc)*(dmiu1*D[f00,miu]+dpz1*D[f00,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjr = 1/I/(w0-k0z*vz-t*wc)*(dmiu0*D[f00,miu]+dpz0*D[f00,pz])/.{t->nr}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};

fs1QLn1 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1Lr)+dy0*(I*k1y*f1Lr)+dth0*(I*S*r*f1Lr)+dmiu0*D[f1Lr,miu]+dpz0*D[f1Lr,pz])/.{t->nr};
fs1QLn2 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx1*(I*k0x*f0Lnjr)+dy1*(I*k0y*f0Lnjr)+dth1*(I*S*nr*f0Lnjr)+dmiu1*D[f0Lnjr,miu]+dpz1*D[f0Lnjr,pz])/.{t->r};

fs1QL1nA = fs1QLn1/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
fs1QL2nA = fs1QLn2/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

fs1QL1nB = fs1QL1nA/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,
BesselJ[nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-nr wc)->b0nr};
fs1QL2nB = fs1QL2nA/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,
BesselJ[nr,a0 x]->j0nr0,BesselJ[-1+nr,a0 x]->j0nr1,BesselJ[1+nr,a0 x]->j0nr2,BesselJ[-2+nr,a0 x]->j0nr3,BesselJ[2+nr,a0 x]->j0nr4}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-nr wc)->b0nr};

(*special simplify*)
fs1QL1nC = fs1QL1nB/.{a0x->a0,a0y->0};
fs1QL2nC = fs1QL2nB/.{a0x->a0,a0y->0};

NA1 = fs1QL1nC*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NA2 = fs1QL2nC*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NB1=NA1*Sqrt[Pi]*Exp[u^2+x^2]/.{ph1->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2=NA2*Sqrt[Pi]*Exp[u^2+x^2]/.{ph1->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1=Simplify[NB1/(wp2*yita0/ns^2/vt^2)]
NC2=Simplify[NB2/(wp2*yita0/ns^2/vt^2)]






