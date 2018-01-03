/*********************************************************************
FGMM: a fast(er) and light Gaussian mixture model implementation.
Copyright (C) 2010  Florent D'Halluin
Contact: florent.dhalluin@epfl.ch

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License,
version 3 as published by the Free Software Foundation.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free
Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*********************************************************************/
#include "gaussian.h"
//#include "smat.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>

/* check the inverse covariance computation */ 
/* #define CHECK_INVERSE  */ 


void dump(struct gaussian* g)
{
	int k=0;
	printf("  prior : %f \n",g->prior);
	printf("  mean : ");
	for(k=0;k<g->dim;k++)
		printf("%f  ",g->mean[k]);
	printf("\n");

	printf("  covariance : ");
	/*for(k=0;k<6;k++)
	printf("%f  ",g->covar[k]);*/
	smat_pmatrix(g->covar);
}

void invert_covar(struct gaussian* g)
{
  float det=1.;
  int i=0,j=0;
  float * pichol, * chol;
  if(!smat_cholesky(g->covar,g->covar_cholesky))
    {
      // g->covar is not full ranked .. adding some noise
      smat_add_diagonal(g->covar, 1.);
      // if after that the covariance is still not positive, we are into
      // big big trouble so let's just give up here.. something went horribly 
      // wrong before .. 
      assert(smat_cholesky(g->covar,g->covar_cholesky));
    }
  pichol = g->icovar_cholesky->_;
  chol = g->covar_cholesky->_;

  for(i=0;i<g->dim;i++)
    {
      det *= *chol;
      *pichol = 1./(*chol);

      chol++;
      pichol++;

      for(j=i+1;j<g->dim;j++)
	{
	  *pichol++ = *chol++;
	}
    }

  det = det*det;
	g->nfactor = sqrtf( pow(2.0 * M_PI, g->dim) * det);
  //g->nfactor = sqrtf( pow(2 * M_PI, g->dim) * det);

  if(g->nfactor <= FLT_MIN)
    {
      // almost non invertible gaussian :: lets add some noise
      g->nfactor = FLT_MIN;
      smat_add_diagonal(g->covar, 1.);
      printf("determinant :: %e\n", det);
      invert_covar(g);
      //exit(0);
    }
}

void gaussian_init(struct gaussian * g,int dim)
{
	int i;
	g->dim = dim;
	g->mean = (float *) malloc(dim * sizeof(float));
	g->covar = NULL;
	g->covar_cholesky = NULL;
	g->icovar_cholesky = NULL;
	smat_zero(&(g->icovar_cholesky),dim);
	for(i=0;i<dim;i++)
		g->mean[i] = 0.;
	smat_zero(&(g->covar),dim);
	smat_identity(g->covar); // just in case :) 
	smat_zero(&(g->covar_cholesky),dim);
	invert_covar(g);
}

void gaussian_free(struct gaussian * g)
{
	free(g->mean);
	smat_free(&g->covar);
	smat_free(&g->covar_cholesky);
	smat_free(&g->icovar_cholesky);
}
/*
void init_random(struct gaussian3d* g)
{
	int k=0;
	for(k=0;k<3;k++)
	{
		g->mean[k] = (float)rand() / RAND_MAX;
		g->covar[k] = 1.;
		g->icovar[k] = 1.;
		g->covar[k+3] = 0.;
		g->icovar[k+3] = 0.;
	}
	invert_covar(g);
}*/


void gaussian_draw(struct gaussian * g, float * out)
{
  int i=0;
  float * tvec;
  tvec = (float *) malloc(g->dim * sizeof(float)); // irk, 
  for(;i<g->dim;i++)
    tvec[i] = randn_boxmuller();
  smat_multv_lt(g->covar_cholesky,tvec,out);
  for(i=0;i<g->dim;i++)
    out[i] += g->mean[i];
  free(tvec);
}

void gaussian_get_subgauss(struct gaussian* g, struct gaussian* result,
						   int n_dim, int * dims)
{
	if(result->dim != n_dim)
	{
		gaussian_free(result);
		gaussian_init(result,n_dim);
	}
	smat_get_submatrix(g->covar,result->covar,n_dim,dims);
	int i=0;
	for(;i<n_dim;i++)
		result->mean[i] = g->mean[dims[i]];
	invert_covar(result);
}

void gaussian_update(struct gaussian * g,
					 const float * data, 
					 float lr)
{
	int i=0;
	int j=0;
	int curs=0;
	for(;i<g->dim;i++)
	{
		g->mean[i] += lr*(data[i] - g->mean[i]);
		for(j=i;j<g->dim;j++)
		{
			g->covar->_[curs] += lr*( (data[i]-g->mean[i])*(data[j] - g->mean[j]) - g->covar->_[curs]) ;
			curs++;
		}
	}
}

