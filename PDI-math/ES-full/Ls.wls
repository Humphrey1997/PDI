#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n00*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};

f0Ln = 1/I/(ws-ksz*vz-t*wc)*(dmius*D[f00,miu]+dpzs*D[f00,pz])/.{t->n};
f0LnA = f0Ln/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.
{ksz->nsz,k1z->n1z,k2z->n2z,ksz->nsz}/.{ksp->as*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{ksx->asx*wc/vt,ksy->asy*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
f0LnB = f0LnA/.{BesselJ[n,as x]->jsn0,BesselJ[n-1,as x]->jsn1,BesselJ[n+1,as x]->jsn2}/.
{nsz->cs/vt,n1z->c1/vt,nsz->c0/vt}/.{(-cs u+ws-n wc)->bsn};

NA = f0LnB*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n00*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NB = NA*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1};
NC = Simplify[NB/(2*wp2/ns^2/vt^2)]
ND = NC/.{jsn0->BesselJ[n,as*x],jsn1->BesselJ[n-1,as*x],jsn2->BesselJ[n+1,as*x]}



