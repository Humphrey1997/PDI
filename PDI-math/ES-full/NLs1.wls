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

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,phj->phs,aj->as};

fsLr=1/I/(ws-ksz*vz-t*wc)*(dmius*D[f00,miu]+dpzs*D[f00,pz])/.{t->r}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lpjrc=-1/I/(w0-k0z*vz-t*wc)*(dmiu0c*D[f00,miu]+dpz0c*D[f00,pz])/.{t->rp}/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

f1QLp1 = 1/I/(w1-k1z*vz-p*wc)*1/2*(dx0c*(I*ksx*fsLr)+dy0c*(I*ksy*fsLr)+dth0c*(I*S*r*fsLr)+dmiu0c*D[fsLr,miu]+dpz0c*D[fsLr,pz])/.{t->rp}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};
f1QLp2 = 1/I/(w1-k1z*vz-p*wc)*1/2*(dxs*(-I*k0x*f0Lpjrc)+dys*(-I*k0y*f0Lpjrc)+dths*(-I*S*rp*f0Lpjrc)+dmius*D[f0Lpjrc,miu]+dpzs*D[f0Lpjrc,pz])/.{t->r}/.{x->vp/vt,vz->pz/m}/.{vp->Sqrt[2*miu*wc/m]};

fsNL1n1 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1QLp1)+dy0*(I*k1y*f1QLp1)+dth0*(I*S*p*f1QLp1)+dmiu0*D[f1QLp1,miu]+dpz0*D[f1QLp1,pz])/.{t->np};
fsNL1n2 = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1QLp2)+dy0*(I*k1y*f1QLp2)+dth0*(I*S*p*f1QLp2)+dmiu0*D[f1QLp2,miu]+dpz0*D[f1QLp2,pz])/.{t->np};

fsNLn1A = fsNL1n1/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3,1/(vt^2 x^2)^(5/2)->1/vt^5/x^5}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};
fsNLn2A = fsNL1n2/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3,1/(vt^2 x^2)^(5/2)->1/vt^5/x^5}/.
{k0z->n0z,k1z->n1z,k2z->n2z,ksz->nsz}/.{k0p->a0*wc/vt,k1p->a1*wc/vt,k2p->a2*wc/vt,ksp->as*wc/vt}/.{k0x->a0x*wc/vt,k0y->a0y*wc/vt,k1x->a1x*wc/vt,k1y->a1y*wc/vt,k2y->a2y*wc/vt,k2x->a2x*wc/vt,ksy->asy*wc/vt,ksx->asx*wc/vt};

fsNLn1B = fsNLn1A/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,BesselJ[-3+r,as x]->jsr5,BesselJ[3+r,as x]->jsr6,
BesselJ[rp,a0 x]->j0rp0,BesselJ[-1+rp,a0 x]->j0rp1,BesselJ[1+rp,a0 x]->j0rp2,BesselJ[-2+rp,a0 x]->j0rp3,BesselJ[2+rp,a0 x]->j0rp4,BesselJ[-3+rp,a0 x]->j0rp5,BesselJ[3+rp,a0 x]->j0rp6,
BesselJ[np,a0 x]->j0np0,BesselJ[-1+np,a0 x]->j0np1,BesselJ[1+np,a0 x]->j0np2}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-p wc)->b1p,(-cs u+ws-r wc)->bsr,(-c0 u+w0-rp wc)->b0rp};
fsNLn2B = fsNLn2A/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,BesselJ[-3+r,as x]->jsr5,BesselJ[3+r,as x]->jsr6,
BesselJ[+rp,a0 x]->j0rp0,BesselJ[-1+rp,a0 x]->j0rp1,BesselJ[1+rp,a0 x]->j0rp2,BesselJ[-2+rp,a0 x]->j0rp3,BesselJ[2+rp,a0 x]->j0rp4,BesselJ[-3+rp,a0 x]->j0rp5,BesselJ[3+rp,a0 x]->j0rp6,
BesselJ[+np,a0 x]->j0np0,BesselJ[-1+np,a0 x]->j0np1,BesselJ[1+np,a0 x]->j0np2}/.
{nsz->cs/vt,n1z->c1/vt,n0z->c0/vt}/.{(-cs u+ws-n wc)->bsn,(-c1 u+w1-p wc)->b1p,(-cs u+ws-r wc)->bsr,(-c0 u+w0-rp wc)->b0rp};

(*special simplify*)
fsNLn1C = fsNLn1B/.{a0x->a0,a0y->0};
fsNLn2C = fsNLn2B/.{a0x->a0,a0y->0};

NA1 = fsNLn1C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NA2 = fsNLn2C*m*vt*m*x*vt/wc*vt*2*Pi*BesselJ[n,as*x]*(S*q)/ep0/ks^2/.{ep0->n0*q*q/m/wp2}/.{ks->ns}/.{BesselJ[n,as*x]->jsn0};
NB1=NA1*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1}/.{ph0->yita0*1/2*m*vt^2/q};
NB2=NA2*Sqrt[Pi]*Exp[u^2+x^2]/.{phs->1}/.{ph0->yita0*1/2*m*vt^2/q};
NC1 = Simplify[NB1/(wp2*yita0^2/ns^2/vt^2)]
NC2 = Simplify[NB2/(wp2*yita0^2/ns^2/vt^2)]

