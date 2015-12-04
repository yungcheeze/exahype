#include "EulerFlow3d/dg/BaseFunc.h"

void exahype::dg::BaseFunc(double phi[EXAHYPE_ORDER+1], double phi_xi[EXAHYPE_ORDER+1], double xi, double xin[EXAHYPE_ORDER+1])
{
	int i, j, m; 
	double tmp;    
	// Initialize variables 
	for(i=0;i<EXAHYPE_ORDER+1;i++)
	{
		phi[i]      = 1.0; 
		phi_xi[i]   = 0.0;
	}
	// Lagrange polynomial and its derivative
	for(m=0;m<EXAHYPE_ORDER+1;m++)
	{
		for(j=0;j<EXAHYPE_ORDER+1;j++)
		{
			if(j==m) continue; 
			phi[m] = phi[m]*(xi-xin[j])/(xin[m]-xin[j]); 
		}
		for(i=0;i<EXAHYPE_ORDER+1;i++)
		{
			if(i==m) continue; 
			tmp = 1.0;  
			for(j=0;j<EXAHYPE_ORDER+1;j++)
			{
				if(j==i) continue;
				if(j==m) continue;
				tmp = tmp*(xi-xin[j])/(xin[m]-xin[j]); 
			} 
			phi_xi[m] = phi_xi[m] + tmp/(xin[m]-xin[i]);  
		} 
	}
}
