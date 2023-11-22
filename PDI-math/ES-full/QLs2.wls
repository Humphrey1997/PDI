#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

dx2 = dx/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dy2 = dy/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dmiu2 = dmiu/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dth2 = dth/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dpz2 = dpz/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};

f2Lr = 1/I/(w2-k2z*vz-t*wc)*(dmiu2*D[f00,miu]+dpz2*D[f00,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
f0Lnjrc = -1/I/(w0-k0z*vz-t*wc)*(dmiu0c*D[f00,miu]+dpz0c*D[f00,pz])/.{t->rn}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

fs2QL1n = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0c*(I*k2x*f2Lr)+dy0c*(I*k2y*f2Lr)+dth0c*(I*S*r*f2Lr)+dmiu0c*D[f2Lr,miu]+dpz0c*D[f2Lr,pz])/.{t->rn};
fs2QL2n = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx2*(-I*k0x*f0Lnjrc)+dy2*(-I*k0y*f0Lnjrc)+dth2*(-I*S*rn*f0Lnjrc)+dmiu2*D[f0Lnjrc,miu]+dpz2*D[f0Lnjrc,pz])/.{t->r};

fs2QL1nA = fs2QL1n/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
fs2QL2nA = fs2QL2n/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

fs2QL1nB = fs2QL1nA/.{BesselJ[r,a2 x]->j2r0,BesselJ[-1+r,a2 x]->j2r1,BesselJ[1+r,a2 x]->j2r2,BesselJ[-2+r,a2 x]->j2r3,BesselJ[2+r,a2 x]->j2r4,
BesselJ[+rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c2 u+w2-r wc)->b2r,(-c0 u+w0-(+rn) wc)->b0rn};
fs2QL2nB = fs2QL2nA/.{BesselJ[r,a2 x]->j2r0,BesselJ[-1+r,a2 x]->j2r1,BesselJ[1+r,a2 x]->j2r2,BesselJ[-2+r,a2 x]->j2r3,BesselJ[2+r,a2 x]->j2r4,
BesselJ[+rn,a0 x]->j0rn0,BesselJ[-1+rn,a0 x]->j0rn1,BesselJ[1+rn,a0 x]->j0rn2,BesselJ[-2+rn,a0 x]->j0rn3,BesselJ[2+rn,a0 x]->j0rn4}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c2 u+w2-r wc)->b2r,(-c0 u+w0-(+rn) wc)->b0rn};

(*special simplify*)
fs2QL1nC = fs2QL1nB/.{a0x->a0,a0y->0};
fs2QL2nC = fs2QL2nB/.{a0x->a0,a0y->0};

NA1 = fs2QL1nC*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NA2 = fs2QL2nC*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NB1=NA1*Sqrt[Pi]*Exp[u^2+x^2]/.{ph2->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2=NA2*Sqrt[Pi]*Exp[u^2+x^2]/.{ph2->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1=Simplify[NB1/(wp2*yita0/ns^2/vt^2)]
NC2=Simplify[NB2/(wp2*yita0/ns^2/vt^2)]









