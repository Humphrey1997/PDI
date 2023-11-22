#!/usr/bin/env wolframscript
(* ::Package:: *)

f00 = n0*wc/m^2/Pi^(3/2)/vt^3*Exp[-(miu*wc+pz^2/2/m)/T];
fsLr = S*q*f00/I/T/(ws-ksz*vz-r*wc)*(vz*Esz*BesselJ[r,as*x]+vp/2*((Esx-I*Esy)*Exp[I*des]*BesselJ[r-S,as*x]+(Esx+I*Esy)*Exp[-I*des]*BesselJ[r+S,as*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
f0Lpjrc = S*q*f00/(-I)/T/(w0-k0z*vz-(-p+r)*wc)*(vz*E0zc*BesselJ[-p+r,a0*x]+vp/2*((E0xc+I*E0yc)*Exp[-I*de0]*BesselJ[-p+r-S,a0*x]+(E0xc-I*E0yc)*Exp[I*de0]*BesselJ[-p+r+S,a0*x]))/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};

dx0 = dx/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dy0 = dy/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dmiu0 = dmiu/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dth0 = dth/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};
dpz0 = dpz/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejx->E0x,Ejy->E0y,Ejz->E0z,aj->a0,wj->w0,dej->de0};

dx0c = dxc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dy0c = dyc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dmiu0c = dmiuc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dth0c = dthc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};
dpz0c = dpzc/.{kjx->k0x,kjy->k0y,kjz->k0z,kjp->k0p,Ejxc->E0xc,Ejyc->E0yc,Ejzc->E0zc,aj->a0,wj->w0,dej->de0};

dxs = dx/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dys = dy/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dmius = dmiu/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dths = dth/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};
dpzs = dpz/.{kjx->ksx,kjy->ksy,kjz->ksz,kjp->ksp,Ejx->Esx,Ejy->Esy,Ejz->Esz,aj->as,wj->ws,dej->des};

f1QLp1 = 1/I/(w1-k1z*vz-p*wc)*1/2*(dx0c*(I*ksx*fsLr)+dy0c*(I*ksy*fsLr)+dth0c*(I*S*r*fsLr)+dmiu0c*D[fsLr,miu]+dpz0c*D[fsLr,pz])/.{t->-p+r};
f1QLp2 = 1/I/(w1-k1z*vz-p*wc)*1/2*(dxs*(-I*k0x*f0Lpjrc)+dys*(-I*k0y*f0Lpjrc)+dths*(I*S*(p-r)*f0Lpjrc)+dmius*D[f0Lpjrc,miu]+dpzs*D[f0Lpjrc,pz])/.{t->r};
f1QLp = (f1QLp1+f1QLp2)/.{vz->pz/m,x->vp/vt}/.{vp->Sqrt[2*miu*wc/m]};
fsNLn = 1/I/(ws-ksz*vz-n*wc)*1/2*(dx0*(I*k1x*f1QLp)+dy0*(I*k1y*f1QLp)+dth0*(I*S*p*f1QLp)+dmiu0*D[f1QLp,miu]+dpz0*D[f1QLp,pz])/.{t->n-p};

fsNLnA = Simplify[fsNLn/.{(E0x+I*E0y)->e0p*beta*m*w0*vt/q,(E0xc-I*E0yc)->e0pc*beta*m*w0*vt/q,(E0x-I*E0y)->e0s*beta*m*w0*vt/q,(E0xc+I*E0yc)->e0sc*beta*m*w0*vt/q}/.
{E0z->0,E0y->e0y*beta*m*w0*vt/q,E0x->e0x*beta*m*w0*vt/q,E0zc->0,E0yc->e0yc*beta*m*w0*vt/q,E0xc->e0xc*beta*m*w0*vt/q}/.{Exp[I*des]->Eis,Exp[-I*des]->Eisc,Exp[I*de0]->Ei0,Exp[-I*de0]->Ei0c,Exp[I*de1]->Ei1,Exp[-I*de1]->Ei1c}/.
{ksz->0,k1z->0,k0z->0,ksx->asx*wc/vt,k1x->a1x*wc/vt,k0x->a0x*wc/vt,ksy->asy*wc/vt,k1y->a1y*wc/vt,k0y->a0y*wc/vt,ksp->as*wc/vt,k1p->a1*wc/vt,k0p->a0*wc/vt}/.{T->1/2*m*vt^2,pz->m*vz,miu->m*vp^2/2/wc}/.{vp->x*vt,vz->u*vt}/.
{Sqrt[vt^2 x^2]->vt*x,1/Sqrt[vt^2 x^2]->1/vt/x,1/(vt^2 x^2)^(3/2)->1/vt^3/x^3}/.{(ws-n wc)->bsn,(w1-p wc)->b1p,(ws-r wc)->bsr, (w0-(-p+r) wc)->b0rp}];

fsNLnB=fsNLnA/.{BesselJ[r,as x]->jsr0,BesselJ[-1+r,as x]->jsr1,BesselJ[1+r,as x]->jsr2,BesselJ[-2+r,as x]->jsr3,BesselJ[2+r,as x]->jsr4,BesselJ[-3+r,as x]->jsr5,BesselJ[3+r,as x]->jsr6,
BesselJ[r-p,a0 x]->j0rp0,BesselJ[-1+r-p,a0 x]->j0rp1,BesselJ[1+r-p,a0 x]->j0rp2,BesselJ[-2+r-p,a0 x]->j0rp3,BesselJ[2+r-p,a0 x]->j0rp4,BesselJ[-3+r-p,a0 x]->j0rp5,BesselJ[3+r-p,a0 x]->j0rp6,
BesselJ[n-p,a0 x]->j0np0,BesselJ[-1+n-p,a0 x]->j0np1,BesselJ[1+n-p,a0 x]->j0np2};
(*special simplify*)
fsNLnC = fsNLnB/.{a0x->a0,a0y->0,de0->0,Ei0->1,Ei0c->1}

Jxy= Simplify[fsNLnC*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*(vp/2)*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt}/.{Sqrt[vt^2 x^2]->vt*x,vz->vt*u,1/Sqrt[vt^2 x^2]->1/vt/x}];
Jzz= Simplify[fsNLnC*m*m*vp/wc*(I/ep0/ws)*(-2*Pi*S*q)*vz*vt*vt/.{ep0->n0*q*q/m/wp2,vp->x*vt}/.{Sqrt[vt^2 x^2]->vt*x,vz->vt*u,1/Sqrt[vt^2 x^2]->1/vt/x}];
cmx1 = Jxy/(beta^2 E^-x^2 wp2 )*b0rp*bsn*b1p*bsr/.{Esx->1,Esy->0,Esz->0};
cmx2 = Simplify[Integrate[cmx1,{u,-Infinity,Infinity}]]
cmy1 = Jxy/(beta^2 E^-x^2 wp2 )*b0rp*bsn*b1p*bsr/.{Esx->0,Esy->1,Esz->0};
cmy2 = Simplify[Integrate[cmy1,{u,-Infinity,Infinity}]]
czz1 = Jzz/(beta^2 E^-x^2 wp2 )*b0rp*bsn*b1p*bsr/.{Esx->0,Esy->0,Esz->1};
czz2 = Simplify[Integrate[czz1,{u,-Infinity,Infinity}]]






