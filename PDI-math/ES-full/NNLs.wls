#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,phj->ph0,aj->a0};

dx1 = dx/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dy1 = dy/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dmiu1 = dmiu/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dth1 = dth/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};
dpz1 = dpz/.{kjx->k1x,kjy->k1y,kjz->k1z,kjp->k1p,phj->ph1,aj->a1};

f1Lr=1/I/(w1-k1z*vz-t*wc)*(dmiu1*D[f00,miu]+dpz1*D[f00,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
f0Lpjr=1/I/(w0-k0z*vz-t*wc)*(dmiu0*D[f00,miu]+dpz0*D[f00,pz])/.{t->pr}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};

fs1QLp1 = 1/I/(ws-ksz*vz-p*wc)*1/2*(dx0*(I*k1x*f1Lr)+dy0*(I*k1y*f1Lr)+dth0*(I*S*r*f1Lr)+dmiu0*D[f1Lr,miu]+dpz0*D[f1Lr,pz])/.{t->pr}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
fs1QLp2 = 1/I/(ws-ksz*vz-p*wc)*1/2*(dx1*(I*k0x*f0Lpjr)+dy1*(I*k0y*f0Lpjr)+dth1*(I*S*pr*f0Lpjr)+dmiu1*D[f0Lpjr,miu]+dpz1*D[f0Lpjr,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};

f1NLq1 = 1/I/(w1-k1z*vz-q*wc)*1/2*(dx0c*(I*ksx*fs1QLp1)+dy0c*(I*ksy*fs1QLp1)+dth0c*(I*S*p*fs1QLp1)+dmiu0c*D[fs1QLp1,miu]+dpz0c*D[fs1QLp1,pz])/.{t->pq}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
f1NLq2 = 1/I/(w1-k1z*vz-q*wc)*1/2*(dx0c*(I*ksx*fs1QLp2)+dy0c*(I*ksy*fs1QLp2)+dth0c*(I*S*p*fs1QLp2)+dmiu0c*D[fs1QLp2,miu]+dpz0c*D[fs1QLp2,pz])/.{t->pq}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};

fsNNLn1 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1NLq1)+dy0*(I*k1y*f1NLq1)+dth0*(I*S*q*f1NLq1)+dmiu0*D[f1NLq1,miu]+dpz0*D[f1NLq1,pz])/.{t->nq};
fsNNLn2 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1NLq2)+dy0*(I*k1y*f1NLq2)+dth0*(I*S*q*f1NLq2)+dmiu0*D[f1NLq2,miu]+dpz0*D[f1NLq2,pz])/.{t->nq};

fsNNLn1A = fsNNLn1/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3,1/(vt^2 x^2)^(5/2)->1/vt^5/x^5}/.{E^((2 (-(1/2) m u^2 vt^2-1/2 m vt^2 x^2))/(m vt^2))->1}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
fsNNLn2A = fsNNLn2/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3,1/(vt^2 x^2)^(5/2)->1/vt^5/x^5}/.{E^((2 (-(1/2) m u^2 vt^2-1/2 m vt^2 x^2))/(m vt^2))->1}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

fsNNLn1B = fsNNLn1A/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,BesselJ[-3+r,a1 x]->j1r5,BesselJ[3+r,a1 x]->j1r6,BesselJ[-4+r,a1 x]->j1r7,BesselJ[4+r,a1 x]->j1r8,
BesselJ[pr,a0 x]->j0pr0,BesselJ[-1+pr,a0 x]->j0pr1,BesselJ[1+pr,a0 x]->j0pr2,BesselJ[-2+pr,a0 x]->j0pr3,BesselJ[2+pr,a0 x]->j0pr4,BesselJ[-3+pr,a0 x]->j0pr5,BesselJ[3+pr,a0 x]->j0pr6,BesselJ[-4+pr,a0 x]->j0pr7,BesselJ[4+pr,a0 x]->j0pr8,
BesselJ[pq,a0 x]->j0pq0,BesselJ[-1+pq,a0 x]->j0pq1,BesselJ[1+pq,a0 x]->j0pq2,BesselJ[-2+pq,a0 x]->j0pq3,BesselJ[2+pq,a0 x]->j0pq4,
BesselJ[nq,a0 x]->j0nq0,BesselJ[-1+nq,a0 x]->j0nq1,BesselJ[1+nq,a0 x]->j0nq2}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-q wc)->b1q,(-cs u+ws-p wc)->bsp,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-pr wc)->b0pr};
fsNNLn2B = fsNNLn2A/.{BesselJ[r,a1 x]->j1r0,BesselJ[-1+r,a1 x]->j1r1,BesselJ[1+r,a1 x]->j1r2,BesselJ[-2+r,a1 x]->j1r3,BesselJ[2+r,a1 x]->j1r4,BesselJ[-3+r,a1 x]->j1r5,BesselJ[3+r,a1 x]->j1r6,BesselJ[-4+r,a1 x]->j1r7,BesselJ[4+r,a1 x]->j1r8,
BesselJ[pr,a0 x]->j0pr0,BesselJ[-1+pr,a0 x]->j0pr1,BesselJ[1+pr,a0 x]->j0pr2,BesselJ[-2+pr,a0 x]->j0pr3,BesselJ[2+pr,a0 x]->j0pr4,BesselJ[-3+pr,a0 x]->j0pr5,BesselJ[3+pr,a0 x]->j0pr6,BesselJ[-4+pr,a0 x]->j0pr7,BesselJ[4+pr,a0 x]->j0pr8,
BesselJ[pq,a0 x]->j0pq0,BesselJ[-1+pq,a0 x]->j0pq1,BesselJ[1+pq,a0 x]->j0pq2,BesselJ[-2+pq,a0 x]->j0pq3,BesselJ[2+pq,a0 x]->j0pq4,
BesselJ[nq,a0 x]->j0nq0,BesselJ[-1+nq,a0 x]->j0nq1,BesselJ[1+nq,a0 x]->j0nq2}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-q wc)->b1q,(-cs u+ws-p wc)->bsp,(-c1 u+w1-r wc)->b1r,(-c0 u+w0-pr wc)->b0pr};

(*special simplify*)
fsNNLn1C = fsNNLn1B/.{a0x->a0,a0y->0};
fsNNLn2C = fsNNLn2B/.{a0x->a0,a0y->0};

NA1 = fsNNLn1C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NA2 = fsNNLn2C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NB1=NA1*Sqrt[Pi]/.{ph1->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2=NA2*Sqrt[Pi]/.{ph1->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1=Simplify[NB1/(wp2*yita0^3/ns^2/vt^2),TimeConstraint->10]
NC2=Simplify[NB2/(wp2*yita0^3/ns^2/vt^2),TimeConstraint->10]


