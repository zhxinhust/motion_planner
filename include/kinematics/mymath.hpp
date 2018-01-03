#ifndef _MYMATH_
#define _MYMATH_

#include <cstdlib>
#include <ctime>
#include <Eigen/Dense>
#include <kdl/frames.hpp>

//using namespace Eigen;

class Interval
{
public:
    double up;
    double lower;

    Interval();
    Interval(double a, double b);
    ~Interval();
};

double getrandom_interval(Interval interval);
double getrandom(double a, double b);
double absv(double x);
double get_sign(double a);
Eigen::Matrix3d hatm(Eigen::Vector3d a);

bool isInRange(double x, double a, double b);

Eigen::Matrix4d POE(Eigen::Vector4d w, Eigen::Vector4d p, double theta);
//Matrix4d POE(Vector3d w, Vector3d p, double theta);

Eigen::Matrix3d calEulerTranMat(KDL::Rotation R);

void arr2vector(float *arr, Eigen::VectorXd &vector);
void vector2arr(Eigen::VectorXd vector, float *arr);

#endif
