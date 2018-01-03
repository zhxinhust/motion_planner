//
// Created by zhxin on 16-12-19.
//
#pragma once

#ifndef DUAL_ARM_ROBOT_MYGMM_H
#define DUAL_ARM_ROBOT_MYGMM_H

#include "mymath.hpp"

#include "dual_arm_robot.hpp"

#include <stdio.h>
#include "fgmm++.hpp"
#include "ur3_kinematics_kdl.hpp"


#define PI 3.1415926
#define NBSTATES    35
#define DIMENSION   6

using namespace KDL;
using namespace sensor_msgs;
using namespace Eigen;

void normalizedata(float *data_in, float *preprocessing, int *randcol, int dim, float *data_out);
//void getArmFrameFromData(float *data, Frame & frame_left, Frame & frame_right);
//void getArmFrameFromCenter(Frame frame_center, float r, Frame &frame_left, Frame &frame_right);
//void getCenterFrameFromData(float *data, Frame & frame_center);
//inline void copyDataFromArmFrame(Frame frame_arm, float *data);
//void getVelFromCenter(Frame frame_center, float r, Twist vel_center, Twist &vel_left, Twist &vel_right);
void normalizedata(VectorXd vector_in, float *preprocessing, int *randcol, int dim, VectorXd &vector_out);
#endif //DUAL_ARM_ROBOT_MYGMM_H

class Robot_Gmm
{
public:
    Gmm *Gmm_dual;
    Gmm *Gmm_left;
    Gmm *Gmm_right;

    // likelihood阈值
    float thre;
    float thre_left;
    float thre_right;

    Robot_Gmm();
    ~Robot_Gmm();

    bool isValidDualArm(VectorXd vector);
    bool isValidLeftArm(VectorXd vector);
    bool isValidRightArm(VectorXd vector);
    bool isValidAll(VectorXd vector);
private:
};


