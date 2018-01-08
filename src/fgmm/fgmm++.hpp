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
#ifndef _FGMMPP_HPP_
#define _FGMMPP_HPP_
/**
 * ------------
 * FGMM library 
 * ------------
 *
 * a fast(er) and light Gaussian mixture model implementation. 
 *        C++ bindings 
 *
 * Florent D'halluin <florent.dhalluin@epfl.ch> 
 */


#include <cstdlib>
#include <vector>
#include <iostream>
//extern "C" {

#include "fgmm.h"
#include "gaussian.h"
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
//#include "dual_arm_robot.hpp"
//}

/**
 * Gaussian Mixture Model class 
 */ 
class Gmm
{
public :

  /**
   * @param states : fixed number of states in the model 
   * @param dim    : dimensionnality of the space we are working in
   */

  int dim;
  int ninput;
  int nstates;

  Gmm(int states, int dim)
  {
    //c_gmm = (struct gmm *) malloc(sizeof(struct gmm ));
    fgmm_alloc(&c_gmm,states,dim);
    c_reg = NULL;
    this->dim = dim;
    this->ninput = 0;
    this->nstates = states;
  };

  ~Gmm()
  {
    if(c_reg != NULL) 
	fgmm_regression_free(&c_reg);
    fgmm_free(&c_gmm);
  };

  /**
   * call this before any kind of learning .. 
   * set means of gaussians by picking random points in 
   * the dataset, and set the variance using the variance 
   * of the dataset
   *
   * @param len : # of points in the dataset
   * @param data : dim*len array, datapoints. (row order) 
   */
  void init(float * data,int len)
  {
    fgmm_init_random(c_gmm,data,len);
  };

  /**
  * call this before any kind of learning .. 
  * set means of gaussians by picking random points in 
  * the dataset, and set the variance using the variance 
  * of the dataset
  *
  * @param len : # of points in the dataset
  * @param data : dim*len array, datapoints. (row order) 
  */
  void init(float * data,int len, int initType)
  {
	  switch(initType)
	  {
	  case 0: // random
		  fgmm_init_random(c_gmm,data,len);
		  break;
	  case 1: // uniform
		  fgmm_init_uniform(c_gmm,data,len);
		  break;
	  case 2: // kmeans
		  fgmm_init_kmeans(c_gmm,data,len);
		  break;
	  }
  };


  /**
   * Just print the model's parameter on stdout
   */
  void Dump()
  {
    fgmm_dump(this->c_gmm);
  };

  /**
   * Expectation Maximization Algorithm. 
   */
  int Em(float * data,int len, 
	 float epsilon=1e-4, enum COVARIANCE_TYPE covar_t = COVARIANCE_FULL)
  {
    return fgmm_em(c_gmm,data,len,&likelihood,epsilon,covar_t,NULL);
  };


  float Pdf(const float * obs, float * weights=NULL)
  {
    return fgmm_get_pdf(c_gmm,obs,weights);
  };

  /*float Pdf(Eigen::VectorXd obs, float * weights=NULL)
  {
    int dim = obs.rows();
    float *data = new float(dim);
    vectorxd2random(obs, data, dim);
    float result = fgmm_get_pdf(c_gmm,data,weights);
    delete data;
    return result;
  };
*/
  float Pdf(const float * obs, int state)
  {
	  if(state >= c_gmm->nstates) return 0;
	  return gaussian_pdf(&c_gmm->gauss[state], obs);
  };

  /**
   * set Prior probability for the desired state
   */
  void SetPrior(int state, float val)
  {
    fgmm_set_prior(this->c_gmm,state,val);
  };
  
  /**
   * set the mean of the specified state, 
   *
   * @param state : the state index 
   * @param mean : an array of size dim, specify the mean
   */
  void SetMean(int state, float * mean)
  {
    fgmm_set_mean(this->c_gmm,state,mean);
  };

  /**
   * set the covariance of the specified state 
   *
   * @param  state : the state index
   * @param  covar : covariance matrix 
   * @param  AsSymetric : Using symetric matrix 
   *                 order .. dim*(dim+1)/2 
   * 
   *      Symetric matrix form : 
   *         [[ 1 2 3 4 ] 
   *          [ 2 5 6 7 ]
   *          [ 3 6 8 9 ]
   *          [ 4 7 9 10]]  
   *  
   * if not we are using a standart row order . 
   */
 
  void SetCovariance(int state, float * covar, bool AsSymetric=true)
  {
    if(AsSymetric) 
      fgmm_set_covar_smat(this->c_gmm,state,covar);
    else
      fgmm_set_covar(this->c_gmm,state,covar);
  };


  float GetPrior(int state)
  {
    return fgmm_get_prior(this->c_gmm,state);
  };

  void GetMean(int state, float * output)
  {
    float * pMean = fgmm_get_mean(this->c_gmm,state);
    for(int i=0;i<this->c_gmm->dim;i++)
      output[i] = pMean[i];
  }

  void GetCovariance(int state, float * out,bool AsSymetric=false)
  {
    if(!AsSymetric)
      {
	fgmm_get_covar(this->c_gmm,state,out);
      }
    else 
      {
	float * pC = fgmm_get_covar_smat(this->c_gmm,state);
	for(int i=0;i<this->c_gmm->dim*(this->c_gmm->dim+1)/2;
	    i++)
	  out[i] = pC[i];
      }
  }
  /**
   * draw a random sample from the model
   *
   * @param sample : the output sample, must be alloc'd of 
   *                 size dim. 
   */
  void Draw(float * sample) 
  {
    fgmm_draw_sample(this->c_gmm,sample);
  };


  /**
   * Initilization for Gaussian Mixture Regression 
   *
   * @param ninput : the first ninput dimension are the inputs, 
   *                 remaining dimensions are outputs. 
   */
  void InitRegression(int ninput){
    if( c_reg != NULL) 
      fgmm_regression_free(&c_reg);
    this->ninput = ninput;
    fgmm_regression_alloc_simple(&c_reg,c_gmm,ninput);
    fgmm_regression_init(c_reg);
  };

  /**
   * Perform the regression on one input point : 
   *
   * @param input : the input point (array of ninput) 
   * @param output : alloc'd array to store result. 
   * @param covar : eventually store resulting covariance is symetric matrix
   *                order (set SetCovariance ) 
   */
  void DoRegression(const float * input, float * output, float * covar=NULL)
  {
    fgmm_regression(c_reg,input,output,covar);
  };


  /**
   * Conditional sampling from the model : 
   * draw a sample in output subspace given the input point in 
   * input subspace. 
   * you must call InitRegression before. 
   */
   void DoSamplingRegression(const float * input, float * output)
   {
     fgmm_regression_sampling(c_reg,input,output);
   };


  /**
   * Online learning HIGHLY EXPERIMENTAL :: 
   * 
   * @param point : input point 
   * @param wta   : use winner take all update ( only update 
   *                 Most likely gaussian) 
   */
  void Update(const float * point,bool wta=false)
  {
    if(wta) 
      fgmm_update_wta(c_gmm,point);
    else 
      fgmm_update(c_gmm,point);
  };


  int GetLikelyState(const float * point)
  {
    return fgmm_most_likely_state(this->c_gmm,point);
  };

  /**
   * 直接增量学习
   * @param data
   * @param pix
   * @param nbData
   * @param nbData0
   */
    void directIncrementalLearning(std::vector<VectorXd> data, MatrixXf pix, int nbData, int nbData0)
    {
        double loglik_threshold = 0.00000000001;
        int nbStep = 0;

        //    MatrixXf MuVector;
        //    VectorXf PriorMat;

        float *Mu, *Sigma;
        float *Priors;
        Mu = new float [dim];
        Sigma = new float[dim * dim];
        Priors = new float[nstates];

        // 将path路径缓存到此
        float *dataTemp;
        dataTemp = new float [dim];

        VectorXf PriorVector, MuVector;       // 前验概率
        MatrixXf SigmaMat;  // 均值和方差，方差只能保存一个矩阵
        PriorVector.resize(nstates);
        MuVector.resize(dim);
        SigmaMat.resize(dim, dim);

        VectorXf Prior0Vector;       // 老的前验概率
        MatrixXf Mu0Mat, Sigma0Mat;     // 老的均值矩阵与方差矩阵
        Prior0Vector.resize(nstates);
        Mu0Mat.resize(dim, nstates);
        Sigma0Mat.resize(dim, dim);

        std::vector<MatrixXf> Sigma0Ful;    // 保存所有的方差矩阵

        // 将初始值保存下来
        for(int i = 0; i < nstates; i++)
        {
            Prior0Vector[i] = GetPrior(i);
            GetMean(i, Mu);
            GetCovariance(i, Sigma, false);

            for(int kk = 0; kk < dim; kk++)
            {
                Mu0Mat(kk, i) = Mu[kk];
                for(int tt = 0; tt < dim; tt++)
                    Sigma0Mat(kk, tt) = Sigma[kk * dim + tt];
            }
            Sigma0Ful.push_back(Sigma0Mat);
        }

        MatrixXf PixMat;    //
        MatrixXf PxiMat;    // 在每一个高斯中的
        PxiMat.resize(nbData, nstates);
        PixMat.resize(nbData, nstates);

        // 将data中的数据保存到矩阵中，方便计算时使用
        MatrixXf DataMat;
        DataMat.resize(dim, nbData);
        for(int i = 0; i < nbData; i++)
            for(int j = 0; j < dim; j++)
                DataMat(j, i) = (float)data[i][j];

        // 后验概率
        VectorXf E0, E;
        E0.resize(nstates);
        E.resize(nstates);

        // 对后验概率累计得到第i个高斯的值
        for(int i = 0; i < nstates; i++)
          E0[i] = pix.col(i).sum();

        std::cout << pix << std::endl;
        std::cout << E0.transpose() << std::endl;

        // 将先验概率保存到向量中，方便后面计算使用
        PriorVector = Prior0Vector;

        MatrixXf tempVector;
        tempVector.resize(1, nstates);

        MatrixXf covTemp;
        covTemp.resize(dim, dim);

        VectorXf F, logF;
        F.resize(nbData);
        logF.resize(nbData);

        float meanLog, meanLog_old = INFINITY;

        while(1)
        {
            /*************** E-step **************/

            // 计算p(x|i)
            for(int i = 0; i < nstates; i++)
            {
                for(int j = 0; j < nbData; j++)
                {
                    // 将第j组数据保存到数组里
                    for(int k = 0; k < dim; k++)
                        dataTemp[k] = (float)data[j][k];

             //       std::cout << data[j].transpose() << std::endl;
                    // 求第j组数据在第i个高斯分量的似然概率
                    PxiMat(j, i) = Pdf(dataTemp, i);
                }
            }


          //  std::cout<<PxiMat<<std::endl;
            // 计算后验概率p(i|x)
            for(int j = 0; j < nbData; j++)
            {
           //     std::cout << "Prior0:" << Prior0Vector.transpose() << std::endl << "PxiMat"<< PxiMat.row(j) << std::endl;
                tempVector = PriorVector.transpose().cwiseProduct(PxiMat.row(j));
//                tempVector = PriorVector.cwiseProduct(PxiMat.row(j));
                PixMat.row(j) = tempVector / tempVector.sum();
            }

            for(int i = 0; i < nstates; i++)
                E[i] = PixMat.col(i).sum();

            std::cout << E.transpose() << std::endl;

            for(int i = 0; i < nstates; i++)
            {
                // 更新先验概率
                PriorVector[i] = (E0[i] + E[i]) / (nbData + nbData0);
                // 更新均值
                MuVector = (Mu0Mat.col(i)*E0(i) + DataMat * PixMat.col(i)) / (E0[i] + E[i]);

                std::cout << "MuVector: " <<MuVector.transpose() << std::endl;

                covTemp.setZero();
                for(int j = 0; j < nbData; j++)
                    covTemp = covTemp + (DataMat.col(j) - MuVector) * (DataMat.col(j) - MuVector).transpose() * PixMat(j, i);

                SigmaMat = ((Sigma0Ful[i] + (Mu0Mat.col(i) - MuVector) * (Mu0Mat.col(i) - MuVector).transpose()) * E0[i] + covTemp) / (E0[i] + E[i]);

                std::cout << "SigmaMat:" <<SigmaMat << std::endl;
                for(int j = 0; j < dim; j++)
                {
                    for(int k = 0; k < dim; k++)
                        Sigma[j * dim + k] = SigmaMat(j, k);

                    Mu[j] = MuVector[j];
                }
                // 设置新的值
                SetCovariance(i, Mu, false);
                SetMean(i, Mu);
                SetPrior(i, PriorVector[i]);
            }

            // 计算p(x|i)
            for(int i = 0; i < nstates; i++)
            {
                for(int j = 0; j < nbData; j++)
                {
                    // 将第j组数据保存到数组里
                    for(int k = 0; k < dim; k++)
                        dataTemp[k] = (float)data[j][k];

                    // 求第j组数据在第i个高斯分量的似然概率
                    PxiMat(j, i) = Pdf(dataTemp, i);
                }
            }

            F = PxiMat * PriorVector;
            logF = Eigen::log(F.array());
            meanLog = logF.mean();

            if(fabs(meanLog / meanLog_old - 1.0) < loglik_threshold)
                break;

            meanLog_old = meanLog;
            nbStep++;
        }

        delete[](Mu);
        delete[](Sigma);
        delete[](Priors);
        delete[](dataTemp);
    };

    struct gmm * c_gmm;
private :

  struct fgmm_reg * c_reg;
  float likelihood;

};

#endif // _FGMMPP_HPP_
