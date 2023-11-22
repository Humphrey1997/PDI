#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n00*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

f0Ln = 1/I/(w0-k0z*vz-t*wc)*(dmiu0*D[f00,miu]+dpz0*D[f00,pz])/.{t->n};
f0LnA = f0Ln/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
f0LnB = f0LnA/.{BesselJ[n,a0 x]->j0n0,BesselJ[n-1,a0 x]->j0n1,BesselJ[n+1,a0 x]->j0n2}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-c0 u+w0-n wc)->b0n};

NA = f0LnB*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a0*x]*(S*q)/ep0/k0^2/.{ep0->n00*q*q/m/wp2}/.{k0->n0}/.{BesselJ[n,a0*x]->j0n0};
NB = NA*Sqrt[Pi]*Exp[u^2+x^2]/.{ph0->1};
NC = Simplify[NB/(2*wp2/n0^2/vt^2)]
ND = NC/.{j0n0->BesselJ[n,a0*x],j0n1->BesselJ[n-1,a0*x],j0n2->BesselJ[n+1,a0*x]}
