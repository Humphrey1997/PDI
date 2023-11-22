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

dx2 = dx/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dy2 = dy/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dmiu2 = dmiu/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dth2 = dth/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};
dpz2 = dpz/.{kjx->k2x,kjy->k2y,kjz->k2z,kjp->k2p,phj->ph2,aj->a2};

f2Lr = 1/I/(w2-k2z*vz-t*wc)*(dmiu2*D[f00,miu]+dpz2*D[f00,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
f0Lpjrc = -1/I/(w0-k0z*vz-t*wc)*(dmiu0c*D[f00,miu]+dpz0c*D[f00,pz])/.{t->rp}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

fs2QLp1 = 1/I/(ws-ksz*vz-p*wc)*1/2*(dx0c*(I*k2x*f2Lr)+dy0c*(I*k2y*f2Lr)+dth0c*(I*S*r*f2Lr)+dmiu0c*D[f2Lr,miu]+dpz0c*D[f2Lr,pz])/.{t->rp}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
fs2QLp2 = 1/I/(ws-ksz*vz-p*wc)*1/2*(dx2*(-I*k0x*f0Lpjrc)+dy2*(-I*k0y*f0Lpjrc)+dth2*(-I*S*rp*f0Lpjrc)+dmiu2*D[f0Lpjrc,miu]+dpz2*D[f0Lpjrc,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};

f2NLn1 = 1/I/(w2-k2z*vz-n*wc)*1/2*(dx0*(I*ksx*fs2QLp1)+dy0*(I*ksy*fs2QLp1)+dth0*(I*S*p*fs2QLp1)+dmiu0*D[fs2QLp1,miu]+dpz0*D[fs2QLp1,pz])/.{t->np};
f2NLn2 = 1/I/(w2-k2z*vz-n*wc)*1/2*(dx0*(I*ksx*fs2QLp2)+dy0*(I*ksy*fs2QLp2)+dth0*(I*S*p*fs2QLp2)+dmiu0*D[fs2QLp2,miu]+dpz0*D[fs2QLp2,pz])/.{t->np};

f2NLn1A = f2NLn1/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3,1/(vt^2 x^2)^(5/2)->1/vt^5/x^5}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
f2NLn2A = f2NLn2/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3,1/(vt^2 x^2)^(5/2)->1/vt^5/x^5}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

f2NLn1B = f2NLn1A/.{BesselJ[r,a2 x]->j2r0,BesselJ[-1+r,a2 x]->j2r1,BesselJ[1+r,a2 x]->j2r2,BesselJ[-2+r,a2 x]->j2r3,BesselJ[2+r,a2 x]->j2r4,BesselJ[-3+r,a2 x]->j2r5,BesselJ[3+r,a2 x]->j2r6,
BesselJ[rp,a0 x]->j0rp0,BesselJ[-1+rp,a0 x]->j0rp1,BesselJ[1+rp,a0 x]->j0rp2,BesselJ[-2+rp,a0 x]->j0rp3,BesselJ[2+rp,a0 x]->j0rp4,BesselJ[-3+rp,a0 x]->j0rp5,BesselJ[3+rp,a0 x]->j0rp6,
BesselJ[np,a0 x]->j0np0,BesselJ[-1+np,a0 x]->j0np1,BesselJ[1+np,a0 x]->j0np2}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-c2 u+w2-n wc)->b2n,(-cs u+ws-p wc)->bsp,(-c2 u+w2-r wc)->b2r,(-c0 u+w0-rp wc)->b0rp};
f2NLn2B = f2NLn2A/.{BesselJ[r,a2 x]->j2r0,BesselJ[-1+r,a2 x]->j2r1,BesselJ[1+r,a2 x]->j2r2,BesselJ[-2+r,a2 x]->j2r3,BesselJ[2+r,a2 x]->j2r4,BesselJ[-3+r,a2 x]->j2r5,BesselJ[3+r,a2 x]->j2r6,
BesselJ[rp,a0 x]->j0rp0,BesselJ[-1+rp,a0 x]->j0rp1,BesselJ[1+rp,a0 x]->j0rp2,BesselJ[-2+rp,a0 x]->j0rp3,BesselJ[2+rp,a0 x]->j0rp4,BesselJ[-3+rp,a0 x]->j0rp5,BesselJ[3+rp,a0 x]->j0rp6,
BesselJ[np,a0 x]->j0np0,BesselJ[-1+np,a0 x]->j0np1,BesselJ[1+np,a0 x]->j0np2}/.
{nsz->cs/vt,n2z->c2/vt,n0z->c0/vt}/.{(-c2 u+w2-n wc)->b2n,(-cs u+ws-p wc)->bsp,(-c2 u+w2-r wc)->b2r,(-c0 u+w0-rp wc)->b0rp};

(*special simplify*)
f2NLn1C = f2NLn1B/.{a0x->a0,a0y->0};
f2NLn2C = f2NLn2B/.{a0x->a0,a0y->0};

NA1 = f2NLn1C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a2*x]*(S*q)/ep0/k2^2/.{ep0->n0*q*q/m/wp2}/.{k2->n2}/.{BesselJ[n,a2*x]->j2n0};
NA2 = f2NLn2C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,a2*x]*(S*q)/ep0/k2^2/.{ep0->n0*q*q/m/wp2}/.{k2->n2}/.{BesselJ[n,a2*x]->j2n0};
NB1=NA1*Sqrt[Pi]*Exp[u^2+x^2]/.{ph2->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2=NA2*Sqrt[Pi]*Exp[u^2+x^2]/.{ph2->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1 = Simplify[NB1/(wp2*yita0^2/n2^2/vt^2)]
NC2 = Simplify[NB2/(wp2*yita0^2/n2^2/vt^2)]










