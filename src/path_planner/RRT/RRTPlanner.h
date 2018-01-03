//
// Created by zhxin on 17-8-16.
//

#ifndef DUAL_ARM_ROBOT_RRTPLANNER_H
#define DUAL_ARM_ROBOT_RRTPLANNER_H

#define MAXEXTENDTIMES 100000

#include "Planner.h"

enum ExtendTree_result
{
    extend_start_fail,    // 扩展树失败
    extend_end_fail,
    extend_success, // 成功扩展树
    extend_finish
};


class RRTPlanner:public Planner
{
public:

    RRTPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit);
   // RRTPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit, ros::Publisher pub);
    ~RRTPlanner();

    virtual void Init(VectorXd str, VectorXd goal);

    PathSearch_result plan();

    void simplifyPath();

    std::vector<Path_Node> tree_start;

    std::vector<Path_Node> tree_goal;

protected:

    double step;

    VectorXd randomRangeLow;   // 采样范围

    VectorXd randomRangeUp;

    Connection_result connectCSpaceRRT(std::vector<Path_Node> &tree, VectorXd goal, VectorXd &vector_out);

    Connection_result connectCSpaceRRT(VectorXd str, VectorXd goal);

    void addTreeNode(std::vector<Path_Node> &trees, int parent, VectorXd v);

    ExtendTree_result extendCSpaceTreesRRT(std::vector<Path_Node> &str_tree, std::vector<Path_Node> &goal_tree);

    ExtendTree_result extendCSpaceTreesRRTSingle(std::vector<Path_Node> &str_tree, VectorXd goal);

    int findNearestNode(VectorXd v, std::vector<Path_Node> tree);

    double rearrangePath();
};

#endif //DUAL_ARM_ROBOT_RRTPLANNER_H
