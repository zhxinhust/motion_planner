//
// Created by zhxin on 17-8-23.
//

#ifndef DUAL_ARM_ROBOT_GMMGUIDEDPLANNER_H
#define DUAL_ARM_ROBOT_GMMGUIDEDPLANNER_H


#include "fgmm++.hpp"

#include "../RRT/RRTPlanner.h"

//using namespace Eigen;

class DisClass
{
public:
    DisClass()
    {
        index = 0;
        dis = 0;
    }
    ~DisClass() = default;

    int index;
    double dis;
};

class myGmm
{
public:
    Gmm *CSpaceGmm;

    double thr;

    int findNearestGaussian(VectorXd p);
    VectorXd angleDis(VectorXd p);

    double dis(VectorXd p);

    VectorXd MahalanobisDis(VectorXd p);

    bool checkValidationRight(VectorXd p);

    myGmm();
    ~myGmm(){};
};


class GMMGuidedPlanner:public RRTPlanner
{
public:
    myGmm GMM;

//    GMMGuidedPlanner();
    GMMGuidedPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit):RRTPlanner(a, robMoveit)
    {
        vexNum = 0;
        edgeNum = 0;
    }

    ~GMMGuidedPlanner(){}

    void Init(VectorXd str, VectorXd goal);

    int findNearestGaussian(VectorXd p);

    void constructRoadMap();

    std::vector<DisClass> findSortedNeighbors(int centerIndex);

    std::vector<DisClass> findSortedNeighborsProb(int centerIndex);

    Connection_result connectCenter(int i, int j);

    VectorXi findPathDigkstra(int strIndex);

    PathSearch_result plan();

    PathSearch_result singleRRTSearch(VectorXd p, int centerIndex, std::vector<Path_Node> &tree);

    PathSearch_result planningInGMM(std::vector<int> centers, std::vector<VectorXd> &path);

protected:
    double calCenterDis(int i, int j);
    MatrixXi E;
    MatrixXf disMat;    // 节点之间的距离矩阵
    int vexNum;     // 节点数量
    int edgeNum;    // 边的数量
    VectorXd getCenter(int i);

    std::vector<VectorXd> planEndPath(VectorXd str, VectorXd goal);
    ExtendTree_result extendCSpaceTreesRRTInGMM(std::vector<Path_Node> &str_tree
            , std::vector<Path_Node> &goal_tree
            , std::vector<int> centers);
};

#endif //DUAL_ARM_ROBOT_GMMGUIDEDPLANNER_H
