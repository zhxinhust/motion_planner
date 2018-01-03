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
#include "fgmm.h"
#include "gaussian.h"

void fgmm_update(struct gmm * gmm,const float * data_point)
{
  int state_i=0;
  float * pxi;
  float nf = 0;
  float eta = .05;
  float weighted_lr;
   
  pxi = (float *) malloc( sizeof(float) * gmm->nstates);
  for(;state_i<gmm->nstates;state_i++)
    {
      pxi[state_i] = gaussian_pdf(&gmm->gauss[state_i], data_point);
      nf += pxi[state_i];
    }
 
  for(state_i=0;state_i<gmm->nstates;state_i++)
    {
      weighted_lr = eta * pxi[state_i]/nf;
      gaussian_update(&gmm->gauss[state_i],data_point,weighted_lr);
      invert_covar(&gmm->gauss[state_i]);
    }
	free(pxi);	
}


/* winner take all approach */

void fgmm_update_wta(struct gmm * gmm,const float * data_point)
{
  int state_i=0;
  int the_state = 0;
  float eta = .05;
  float max_p = 0.;
  float w =  0.;
  for(;state_i<gmm->nstates;state_i++)
    {
      w = gaussian_pdf(&gmm->gauss[state_i], data_point);
      if(w>max_p) 
	{
	  max_p = w;
	  the_state = state_i;
	}
    }
  gaussian_update(&gmm->gauss[the_state],data_point,eta);
  invert_covar(&gmm->gauss[the_state]);
}
  
  
