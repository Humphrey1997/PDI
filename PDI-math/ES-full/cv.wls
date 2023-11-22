#!/usr/bin/env wolframscript
(* ::Package:: *)

S = 1;

kj = {kjx,kjy,kjz};

Mj = {{I*wj/(wj^2-wc^2),-wc/(wj^2-wc^2),0},{wc/(wj^2-wc^2),I*wj/(wj^2-wc^2),0},{0,0,I/wj}};
Ej = -I*kj*phj;
VjL = -S*q/m*Dot[Mj,Ej];

Mjc = {{-I*wj/(wj^2-wc^2),-wc/(wj^2-wc^2),0},{wc/(wj^2-wc^2),-I*wj/(wj^2-wc^2),0},{0,0,-I/wj}};
Ejc = I*kj*phj;
VjLc = -S*q/m*Dot[Mjc,Ejc];



ks = kj/.{kjx->nsx,kjy->nsy,kjz->nsz};
k0 = kj/.{kjx->n0x,kjy->n0y,kjz->n0z};
k1 = kj/.{kjx->n1x,kjy->n1y,kjz->n1z};

V0L = VjL/.{wj->w0,kjx->n0x,kjy->n0y,kjz->n0z,phj->ph0};
VsL = VjL/.{wj->ws,kjx->nsx,kjy->nsy,kjz->nsz,phj->phs};
V0Lc = VjLc/.{wj->w0,kjx->n0x,kjy->n0y,kjz->n0z,phj->ph0};
VsLc = VjLc/.{wj->ws,kjx->nsx,kjy->nsy,kjz->nsz,phj->phs};

M1 = Mj/.{wj->w1};
Ms = Mj/.{wj->ws};

vkv = Simplify[Simplify[Dot[V0L,ks]*VsLc]- Simplify[Dot[VsLc,k0]*V0L]];
V1QL = I/2*Dot[M1,vkv];
V1QLA = Simplify[V1QL/.{ph0->yita0*1/2*m*vt^2/q}/.{phs->1}];
n1A = Simplify[Dot[k1,V1QLA]/w1]*n00;
n1B = n1A*(S*q)/ep0/n1^2/.{ep0->n00*q*q/m/wp2};
n1C = Simplify[n1B/(wp2*yita0*vt^2/n1^2/4)]

