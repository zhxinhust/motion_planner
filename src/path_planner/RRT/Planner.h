//
// Created by zhxin on 17-8-12.
//

#ifndef DUAL_ARM_ROBOT_PLANNER_H
#define DUAL_ARM_ROBOT_PLANNER_H

#include <moveit/move_group_interface/move_group.h>
#include <moveit/planning_scene_interface/planning_scene_interface.h>

#include <moveit/robot_model_loader/robot_model_loader.h>
#include <moveit/robot_model/robot_model.h>
#include <moveit/robot_state/robot_state.h>

#include <moveit/planning_scene/planning_scene.h>

#include <moveit_msgs/DisplayRobotState.h>
#include <moveit_msgs/DisplayTrajectory.h>

#include <moveit_msgs/AttachedCollisionObject.h>
#include <moveit_msgs/CollisionObject.h>
#include <math.h>

#include "mymath.hpp"

#include <boost/date_time.hpp>
#include <trac_ik/trac_ik.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include "dual_arm_robot.hpp"

#include <stdio.h>

#include <tf/transform_broadcaster.h>

using namespace KDL;
using namespace sensor_msgs;
using namespace Eigen;

enum PathSearch_result
{
    pathsearch_finish,
    pathsearch_fail
};

enum Connection_result
{
    connection_success,
    connection_fail,
    connection_mid
};

/**
 * 路径节点类
 */
class Path_Node
{
public:
    VectorXd vector;    // 节点位置向量
    int parent;         // 父节点
    double dis;         // 从根节点到此节点的距离
    std::vector<int> child;//子节点

    Path_Node();
    Path_Node(int dim);
    Path_Node(VectorXd vector,int parent, double dis);
    ~Path_Node();

    void setNode(VectorXd vector,int parent, double dis);

    Path_Node(VectorXd vector_, int parent_, double dis_, std::vector<int> child_);

};

class  Feasible_Path
{
public:
    std::vector<VectorXd> feasiblepath;
    double total_cost;

    Feasible_Path();
    ~Feasible_Path();
    Feasible_Path(std::vector<VectorXd> feasiblepath_, double total_cost_);

    void saveFeasiblePath(std::vector<VectorXd> feasiblepath_, double total_cost_);
};

/**
 * 碰撞检测所需要的类
 */
class Robot_moveit
{
public:
    planning_scene::PlanningScene *planning_scene;

    collision_detection::CollisionRequest *collision_request_left;
    collision_detection::CollisionResult collision_result_left;

    collision_detection::CollisionRequest *collision_request_right;
    collision_detection::CollisionResult collision_result_right;

    const robot_state::JointModelGroup *joint_model_group_left;
    const robot_state::JointModelGroup *joint_model_group_right;

    moveit::planning_interface::MoveGroup *move_group_left;
    moveit::planning_interface::MoveGroup *move_group_right;

    //robot_state::RobotState &check_collision_state;

    Robot_moveit();
    ~Robot_moveit();
};

/**
 * 路径规划器父类
 */
class Planner
{
public:
  //  Planner();
    Planner(robot_state::RobotState &a, Robot_moveit &robMoveit);

    //Planner(robot_state::RobotState &a, Robot_moveit &robMoveit, ros::Publisher &pub);
    ~Planner();

  //  void PlannerInitialize();

    robot_state::RobotState &check_collision_state;

  //  robot_state::RobotState& check_collision_state;

    bool collisionCheckRight(VectorXd p);

    bool collisionCheckLeft(VectorXd p);

    virtual PathSearch_result plan() = 0;

    void FreeMemory();

    void publishJointState(VectorXd jntLeft, VectorXd jntRight);

    std::vector<VectorXd> path; // 规划得到的path路径

    std::vector<VectorXd> optimalPath;//最优路径

protected:

  //  ros::Publisher &pub;

  //  Connection_result connection(VectorXd str, VectorXd goal);  // 尝试连接str和goal

    VectorXd _str;   // 起始路径

    VectorXd _goal;  // 目标位置

    Robot_moveit &robot_moveit;

  //  ros::Publisher *pub;

  //  ros::NodeHandle node_handle;

};


#endif //DUAL_ARM_ROBOT_PLANNER_H
