//
// Created by zhxin on 17-8-17.
//

#include "RRTPlanner.h"

#ifndef DUAL_ARM_ROBOT_RRTSTARTPLANNER_H
#define DUAL_ARM_ROBOT_RRTSTARTPLANNER_H

class RRTStarPlanner:public RRTPlanner
{
public:
    RRTStarPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit);

    ~RRTStarPlanner();

    PathSearch_result plan();

protected:

    double radius;

    ExtendTree_result extendCSpaceTreesRRTStar(std::vector<Path_Node> &str_tree, std::vector<Path_Node> &goal_tree);

    Connection_result connectCSpaceRRTStar(std::vector<Path_Node> &tree, VectorXd goal, VectorXd &vector_out, int nodeIndex);

    void rewirePathNodes(VectorXd vector_new, std::vector<Path_Node> trees, std::vector<int> &neighbors_others,
                         std::vector<int> &rewire_neighbors, int vector_new_index);

    int chooseBestParent(VectorXd vector_new, std::vector<Path_Node> trees, std::vector<int> &neighbors_others);

    void modifyChildCost(std::vector<Path_Node> trees, int nodeIndex, const int rewireNeighborIndex);
};

#endif //DUAL_ARM_ROBOT_RRTSTARTPLANNER_H
