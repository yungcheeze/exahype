#include "SRHDApplication.h"
#include <stdio.h>

/**
 * Everything in this C++ file was stupidly translated from Fortran.
 **/

bool cons2prim(double* V, const double *Q) {
	DEFINE_LABELS;
	double gam = 1.0 / eos_gamma * (eos_gamma-1.0);
	
	/*const*/ double tol = 1e-8, rho_floor = 10e-4, p_floor=10e-5;
	
	if(Q[rho] < 0.) {
		V[rho] = rho_floor;
		V[vx] = V[vy] = V[vz] = 0;
		V[p] = p_floor;
		return SUCCESS;
	}
	
	double d=Q[0], sx=Q[1], sy=Q[2], sz=Q[3], e=Q[4]+d;
	
	/* from srhd3dfortran, PDE.f90:  solve for v2 (Del Zanna et al. (2007) A&A, 473, 11-30) */
	double s2 = sx*sx + sy*sy + sz*sz;
	double b2=0, sb2=0, eps=1e-8;
	
	double w=0, x1=0, x2=1-eps;
	bool rtfailed = false;
	double v2 = rtsafe_c2p_rmhd1(x1,x2,tol,gam,d,e,s2,b2,sb2,w,rtfailed);
	if(rtfailed) {
		V[p] = p_floor;
		V[rho] = rho_floor;
		V[vx] = V[vy] = V[vz] = 0;
		return FAILED;
	} else {
		double den = 1.0/(w+b2);
		V[rho] = d*sqrt(1-v2);
		V[vx] = sx*den;
		V[vy] = sy*den;
		V[vz] = sz*den;
		V[p] = max(1.e-15, gam*(w*(1-v2)-V[rho]));
		return SUCCESS;
	}
}

bool prim2cons(const double* V, double *Q) {
	DEFINE_LABELS;
	double v2 = V[vx]*V[vx] + V[vy]*V[vy] + V[vz]*V[vz];
	
	if(v2 > 1.0) {
		printf("Superluminal velocity in prim2cons");
		return FAILED;
	}
	
	double lf =  1.0 / sqrt(1.0 - v2);
	double gamma1 = eos_gamma/(eos_gamma-1.0);
	double w = V[rho] + gamma1*V[p];
	double ww = w*lf*lf;
	
	Q[rho] = V[rho]*lf;
	Q[vx] = ww*V[vx];
	Q[vy] = ww*V[vy];
	Q[vz] = ww*V[vz];
	Q[p] = ww - V[p] - Q[rho];
	
	return SUCCESS;
}


double rtsafe_c2p_rmhd1(double& x1, double& x2, double& xacc, double& gam, double& d, double& e, double& s2, double& b2, double& sb2, double& w, bool& failed) {
	const int maxit = 200;
	double retval;
	double  fl,fh,df,xh,xl,swap,dxold,dx,f,temp;
	
	failed = false;
	func_c2p_rmhd1(x1,fl,df,gam,d,e,s2,b2,sb2,w);
	if(fl == 0)
		return x1;
	func_c2p_rmhd1(x2,fh,df,gam,d,e,s2,b2,sb2,w);
	if(fh == 0)
		return x2;
	if(fl*fh > 0) {
		failed = true;
		return retval;
	}
	if(fl > 0.) {
		xl = x1;
		xh = x2;
	} else {
		xh   = x1;
		xl   = x2;
		swap = fl;
		fl   = fh;
		fh   = swap;
	}
	retval = .5*(x1+x2);
	dxold = abs(x2-x1);
	dx = dxold;
	func_c2p_rmhd1(retval,f,df,gam,d,e,s2,b2,sb2,w);
	for(int j=1; j<=maxit; j++) {
		if(((retval-xh)*df-f)*((retval-xl)*df-f) >= 0. || abs(2.*f)>abs(dxold*df)) {
			dxold = dx;
			dx = 0.5*(xh-xl);
			retval = xl+dx;
			if(xl==retval) return retval;
		} else {
			dxold = dx;
			dx = f/df;
			temp = retval;
			retval = retval-dx;
			if(temp==retval) return retval;
		}
		if(abs(dx) < xacc) return retval;
		func_c2p_rmhd1(retval,f,df,gam,d,e,s2,b2,sb2,w);
		if(f<0.) {
			xl = retval;
			fl = f;
		} else {
			xh = retval;
			fh = f;
		}
	}
	failed = true;
	return retval;
}

void func_c2p_rmhd1(double& x, double& f, double& df, double& gam, double& d, double& e, double& s2, double& b2, double& sb2, double& w) {
  //
  // This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  // and it corresponds to their choice 3 in Section 3.2
  //

  double v2,rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2;
  
  v2 = x;
  rho = d*sqrt(1.-v2);
  
  c3 = 1.-gam*(1.-v2);
  c2 = gam*rho+.5*b2*(1.+v2)-e;
  c0 = -0.5*sb2;
  
  // For every x=v2, we solve a cubic in W of the form: 
  // c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  // W=y of the paper. If sb=0 ( meaning c0 = 0), 
  // w = -c2/c3 > 0 and dw = 0 in the do loop below. 
  // If -c2/c3 < 0 when sb=0, which makes w=0, 
  // this is a signature that something was wrong before.
  
  if(abs (c0) < 1.0-20) {
     w = -c2 / c3;
  } else {
     w = max ( - c2 / c3, pow( -c0 / c3, 1./3.));
     for(int iter=1; iter<=100; iter++) {
        dw = -((c3*w + c2)*w*w + c0)/((3*c3*w + 2*c2)*w);
        if(abs(dw/w)<1.e-10) return;
        w = w + dw;
     }
  }
  
  dc3   = gam;
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2));
  dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2);
  wb    = w + b2;
  vb2   = sb2 / (w*w);
  f     = wb*wb * v2 - ( 2.0 * w + b2) * vb2 - s2;
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2));
}
