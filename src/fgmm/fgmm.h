/*********************************************************************
FGMM: a fast(er) and light Gaussian mixture model implementation.
Copyright (C) 2010  Florent D'Halluin, Basilio Noris
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

/**
   ------------
   FGMM library 
   ------------
   
   a fast(er) and light Gaussian mixture model implementation. 
   
   Florent D'halluin <florent.dhalluin@epfl.ch> 
 */

/** @file fgmm.h 
 *  @brief Awesome fast Gaussian Mixture Model library. 
 */

/**  @mainpage FGMM Library
 * 
 * This library provides a fast(er) and light GMM implementation with the following 
 * nice features : 
 *    - sampling
 *    - EM learning 
 *    - online update rule (experimental) 
 *    - No dependencies ( seriously, at all !!) 
 *    - faster that the previous C++ implementation 
 *    - pure C implementation (with a C++ wrapper) 
 *    - comes with Matlab/Octave wrappers and a Python wrapper. 
 *       
 *  To learn more about C API see fgmm.h 
 *  
 *  Or you might want some nice OO programming, then check out Gmm class
 */

#ifndef _FGMM_H_
#define _FGMM_H_

#ifndef NULL
#define NULL 0
#endif
/**
 * structure holding one state of the  gmm 
 *
 * see gaussian.h for more details .. 
 */

struct gaussian;

/**
 * structure holding the whole Gaussian Mixture model 
 * this must be initialised by fgmm_alloc 
 *
 * short example  : 
 * \code
 *   struct gmm mygmm;
 *   fgmm_alloc(&mygmm,4,1) // 1D 4 states GMM 
 *   int i = 0;
 *   for(i=0;i<4;i++)
 *     { 
 *       fgmm_set_prior(&mygmm,i,.25); 
 *        fgmm_set_mean(&mygmm,i,i*2.5):
 *        fgmm_set_covar(&mygmm,i,.1):
 *      }
 *    float sample; 
 *    fgmm_draw_sample(&mygmm,&sample);
 *    fgmm_free(&mygmm);
 * \endcode
 */
       
struct gmm {
  /**
   * opaque structure representing a gaussian 
   * you can access to parameters using 
   * fgmm_get_* functions 
   */
  struct gaussian * gauss;
  /**
   * number of states 
   */
  int nstates;
  /**
   * dimensionnality of the space 
   */
  int dim;
};

/**
 *  alloc all the memory needed by the model (alloc all gaussians struct and all ) 
 *  
 *  all gaussians are init'd to zero mean, unity covariance 
 *  zero prior proba
 *
 *  @param **gmm : pointer to an gmm struct  
 *  @param nstates : number of states 
 *  @param dim : dimensionnality of input space 
 */
void fgmm_alloc(struct gmm ** gmm, int nstates, int dim);

/**
 * free everything allocated by the above function 
 */
void fgmm_free(struct gmm ** gmm);

/**
 * initialize the model from the data by :
 *  for each gaussian :
 *     - pick one random data point from data 
 *     - set the mean to this point 
 *     - set covariance to data covariance/nstates
 *     - set prior to 1./nstates
 *
 * the model must be first alloc'd (with fgmm_alloc) 
 */
void fgmm_init_random(struct gmm * gmm,
		     const float * data,
		     int data_len);

/**
* initialize the model from the data by segmenting the space
*  in uniform splits and for each gaussian:
*     - pick a data point in the middle of the corresponding split
*     - set the mean to this point 
*     - set covariance to data covariance/nstates
*     - set prior to 1./nstates
*
* the model must be first alloc'd (with fgmm_alloc) 
*/
void fgmm_init_uniform(struct gmm * gmm,
					  const float * data,
					  int data_len);

void fgmm_init_kmeans( struct gmm * gmm,
					  const float * data,
					  int data_len);

/**
 * set the prior of a given state 
 *
 * @param state : index of the state
 * @param prior : the prior value (be careful to keep sum(priors) = 1. )
 */
void fgmm_set_prior(struct gmm *,int state, float prior);

/**
 * set the mean of a given state 
 *
 * @param state : index of the state
 * @param mean : table of dim float 
 */

void fgmm_set_mean(struct gmm *,int state, const float * mean);

float fgmm_get_prior(struct gmm * gmm, int state);

/**
 * get a pointer to the given state mean 
 * 
 * this is a pointer to the actual mean, so take care .. 
 */

float * fgmm_get_mean(struct gmm * gmm,int state);

/**
 * Set the covariance of a given state 
 *
 * @param state : index of the state 
 * @param covar : The covariance matrix 
 * 
 *      Symetric matrix form : 
 *         [[ 1 2 3 4 ]
 *          [ 2 5 6 7 ]
 *          [ 3 6 8 9 ]
 *          [ 4 7 9 10]]  
 *
 *  covar must be a float table of lenght (dim*dim+1)/2 
 */

void fgmm_set_covar_smat(struct gmm *,int state, const float * covar);

/**
 * Set the covariance of a given state from a square 
 * matrix (row order : [i*dim + j] =  [i][j] )
 */
void fgmm_set_covar(struct gmm * gmm,int state, 
		    const float * square_covar);

/**
 * get the address of the symetric form of the covar matrix
 *
 * this is the actuall address of covariance matrix, all kind 
 * of trouble can happen to you with this .. 
 */
float * fgmm_get_covar_smat(struct gmm * gmm, int state) ;

/**
 * copy the covariance matrix in a square matrix ( row 
 * order ) 
 *
 * way safer than the smat variant .. 
 * @param square_covar : a dim*dim alloc'd float array
 * @param state : the state index we are interested in 
 */
void fgmm_get_covar(struct gmm * gmm, 
		    int state,
		    float * square_covar); /* -> must be alloc'd */ 

/**
 * print the gmm parameters to screen 
 */

void fgmm_dump(struct gmm * gmm);

/**
 * draw one sample from the gmm 
 *
 * @param out : *alloc'd * table of dim 
 *
 */
void fgmm_draw_sample(struct gmm *, float * out);


enum COVARIANCE_TYPE 
  {
    COVARIANCE_FULL,  // full covariance matrix
    COVARIANCE_DIAG,  // diagonal covariance matrix
    COVARIANCE_SPHERE // "sphered" covariance (ie same value over the diagonal
  };

//#define loglikelihood_eps 1e-4

/**
 * EM algorithm 
 *
 * @param GMM : the initialized gmm 
 * @param data : a dim x data_length float table holding training data (row order) 
 *               the k-th training point is data[dim*k]..data[dim*(k+1)] 
 * @param data_length : number of training points 
 * @param end_loglikelihood : will be set to the final loglikelihood 
 * @param likelihood_epsilon : if the loglikelihood variation is below 
 *        this threshold, stops here
 * @param covar_t : how would you like your covariance today ?? 
 * @param weights : if not NULL, will perform the weighted variant of EM. 
 *                  weight must then be an allocated array of 
 *                  size data_length
 *
 * @return : The number of iterations 
 */

int fgmm_em( struct gmm * GMM,
	     const float * data,
	     int data_length, 
	     float * end_loglikelihood,
	     float likelihood_epsilon,
	     enum COVARIANCE_TYPE covar_t,
	     const float * weights);


static int fgmm_em_simple(struct gmm * GMM, const float * data, int data_length)
{
  float nevermind=0;
  return fgmm_em(GMM,data,data_length,
		 &nevermind, 1e-4, COVARIANCE_FULL, NULL);
};


int fgmm_kmeans( struct gmm * GMM,
				const float * data,
				int data_length,
				float epsilon,
				const float * weights);

/**
 * return likelihood of point
 * if weights != NULL , return normalized weights of each gaussian
 */
float fgmm_get_pdf( struct gmm * gmm,
		    const float * point,
		    float * weights);


/**
 * Structure holding stuffs for the regression 
 * 
 * This must be alloc'd using fgmm_regression_alloc* on a 
 *  non allocated pointer. 
 */

struct fgmm_reg;

/**
 * init a fgmm_reg structure for a regression 
 * were the input_len first input dimenstion are 
 * the inputs and the remaining dimensions the output
 */
void fgmm_regression_alloc_simple(struct fgmm_reg ** regression,
				 struct gmm * gmm,
				 int input_len);

/**
 * init a fgmm_reg structure 
 * @param gmm : the model to perform regression on 
 * @param input_len is the number of input dimensions
 * @param input_dim are the indexes of the input dimensions
 * @param output_len is the number of outputs 
 * @param output_dim is their indexes .. 
 */ 
void fgmm_regression_alloc(struct fgmm_reg ** regression,
			  struct gmm * gmm,
			  int input_len, int * input_dim,
			  int output_len, int * output_dim);

/**
 * free all the memory allocated for the regression structure 
 */
void fgmm_regression_free(struct fgmm_reg ** regression);

/**
 * does all intermediate computation, call this 
 * between alloc and regression. Call this also if 
 * the gmm relative to the regression has changed
 *
 * this caches the inverse covariance matrices for instance
 */
void fgmm_regression_init(struct fgmm_reg * reg);

/**
 * does the regression 
 */
void fgmm_regression(struct fgmm_reg * reg, const float * inputs, 
		     float * outputs, float * covar);


/**
 * Conditional sampling
 *
 * draw a sample in the output subspace, given the input point 
 * out ~ p(x | input) 
 */
void fgmm_regression_sampling(struct fgmm_reg * reg,const float * inputs,
			      float * output);


/**
 * return index of the most likely state, given an observation
 */
int fgmm_most_likely_state(struct gmm *, const float * obs);

/**
 * incremental updates, update the model with a new datapoint
 * 
 * Highly experimental .. 
 */
void fgmm_update(struct gmm * gmm, const float * data_point);

/**
 * on-line update, with winner take all (only update the most 
 * likely gaussian ) 
 */
void fgmm_update_wta(struct gmm * gmm, const float * data_point);

#endif // _FGMM_H_
