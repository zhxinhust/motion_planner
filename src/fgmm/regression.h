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
#ifndef _REGRESSION_H_
#define _REGRESSION_H_
#include "gaussian.h"
#include "fgmm.h"

struct gaussian_reg {
  struct gaussian * gauss;
  struct gaussian * subgauss; // input subgaussian Used to compute the weight of this
  int * input_dim;
  int * output_dim;
  int input_len;
  int output_len;
  float * reg_matrix; // store in->out A matrix 
};


struct fgmm_reg {
  struct gmm * model;
  int * input_dim;
  int * output_dim;
  int input_len;
  int output_len;
  struct gaussian_reg * subgauss;
};

#endif // _REGRESSION_H_
