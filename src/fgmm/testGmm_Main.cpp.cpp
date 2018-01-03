//
// Created by zhaoxin on 17-9-11.
//

#include "../path_planning/GMMPlanner/GMMGuidedPlanner.hpp"
#include "../GMM/GMM2.h"

int main()
{
    Gmm GMM(3, 2);

    for(int i = 0; i < 3; i++)
    {
        GMM.SetPrior(i, prior[i]);
        GMM.SetMean(i, Mu[i]);
        GMM.SetCovariance(i, sigma[i], false);
    }

    std::vector<VectorXd> dataVector;
    VectorXd dataTemp;
    dataTemp.resize(2);

    for(int i = 0; i < nbdata; i++)
    {
        dataTemp[0] = inData[i][0];
        dataTemp[1] = inData[i][1];
        dataVector.push_back(dataTemp);
    }

    MatrixXf pix;
    pix.resize(nbdata, 3);

    for(int i = 0; i < nbdata; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            pix(i, j) = Pix[i][j];
        }
    }

    GMM.directIncrementalLearning(dataVector, pix, nbdata, nbdata);

    return 0;
}