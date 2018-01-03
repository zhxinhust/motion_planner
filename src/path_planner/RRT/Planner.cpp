//
// Created by zhxin on 17-8-12.
//

#include "Planner.h"

Path_Node::Path_Node(){}

Path_Node::Path_Node(int dim)
{
    vector.resize(dim);
    dis = 0;
    parent = 0;
    child.clear();

}

Path_Node::Path_Node(VectorXd vector,int parent, double dis)
{
    this->parent = parent;
    this->dis = dis;
    this->vector = vector;
    child.clear();
}

Path_Node::Path_Node(VectorXd vector_, int parent_, double dis_, std::vector<int> child_)
{
    this->parent = parent_;
    this->dis = dis_;
    this->vector = vector_;
    this->child = child_;
}

Path_Node::~Path_Node(){}


/*double Path_Node::distance_r_theta_phi(VectorXd v)
{
    double dist;
    Frame frame_left, frame_right;
    getArmFrameFromVector(v, frame_left, frame_right);
    Vector temp = frame_left.p - frame_right.p;

    dist = temp.Norm();
    return dist;
}*/

void Path_Node::setNode(VectorXd vector, int parent, double dis)
{
    this->vector = vector;
    this->dis = dis;
    this->parent = parent;
    child.clear();
}

Feasible_Path::Feasible_Path() {}
Feasible_Path::~Feasible_Path() {}

Feasible_Path::Feasible_Path(std::vector<VectorXd> feasiblepath_, double total_cost_)
{
    this->feasiblepath = feasiblepath_;
    this->total_cost = total_cost_;
}

void Feasible_Path::saveFeasiblePath(std::vector<VectorXd> feasiblepath_, double total_cost_)
{
    feasiblepath.clear();
    for(int i=0; i<feasiblepath_.size(); i++)
    {
        this->feasiblepath.push_back(feasiblepath_[i]);
    }
//    std::vector<VectorXd>::iterator it;
//    int length = (int)feasiblepath_.size();
//    for(int i=0; i<length; i++)
//    {
//        it = feasiblepath_.begin();
//        feasiblepath.insert(it, feasiblepath_.at(i));
//    }
//    feasiblepath.assign(feasiblepath_.begin(), feasiblepath_.end());

//    this->feasiblepath = feasiblepath_;
    this->total_cost = total_cost_;
}

/**
 * 构造函数
 */
Robot_moveit::Robot_moveit()
{
    static const std::string PLANNING_GROUP_LEFT = "left_arm";
    static const std::string PLANNING_GROUP_RIGHT = "right_arm";

    move_group_left = new moveit::planning_interface::MoveGroup(PLANNING_GROUP_LEFT);
    move_group_right = new moveit::planning_interface::MoveGroup(PLANNING_GROUP_RIGHT);

    joint_model_group_left = move_group_left->getCurrentState()->getJointModelGroup(PLANNING_GROUP_LEFT);
    joint_model_group_right = move_group_right->getCurrentState()->getJointModelGroup(PLANNING_GROUP_RIGHT);

    robot_model_loader::RobotModelLoader robot_model_loader("robot_description");
    robot_model::RobotModelPtr kinematic_model = robot_model_loader.getModel();

    planning_scene = new planning_scene::PlanningScene(kinematic_model);

    collision_request_left = new collision_detection::CollisionRequest();
    collision_request_right = new collision_detection::CollisionRequest();

    collision_request_right->group_name = "right_arm";
    collision_request_left->group_name = "left_arm";
}

Robot_moveit::~Robot_moveit(){}

Planner::Planner(robot_state::RobotState &a, Robot_moveit &robMoveit):check_collision_state(a),robot_moveit(robMoveit)
{

}

//Planner::Planner(robot_state::RobotState &a, Robot_moveit &robMoveit, ros::Publisher &publisher):check_collision_state(a),robot_moveit(robMoveit), pub(publisher) {}

/*Planner::Planner()
{
  //  *pub = node_handle.advertise<sensor_msgs::JointState>("/joint_states", 10);         // 初始化发布器，发布关节位置
    robot_moveit = new Robot_moveit();  // 碰撞检测所需要用到的类

    robot_moveit->planning_scene->getCurrentStateNonConst();

    check_collision_state = robot_moveit->planning_scene->getCurrentStateNonConst();   // 获取碰撞检测所需的变量


 //   node_handle.param("timeout", timeout, 0.005);
 //   node_handle.param("urdf_param", urdf_param, std::string("/robot_description"));
}*/

Planner::~Planner()
{

}

bool Planner::collisionCheckRight(VectorXd p)
{
    // 将要检测的角度位置更新到模型中去
    std::vector<double> joint_pos(6, 0);

    for(unsigned i = 0; i < 6; i++)
        joint_pos[i] = p[i];

    check_collision_state.setJointGroupPositions(robot_moveit.joint_model_group_right, joint_pos);

    // 清除原来计算结果
    robot_moveit.collision_result_right.clear();

    // 碰撞检测
    robot_moveit.planning_scene->checkCollision(*robot_moveit.collision_request_right,
                                                 robot_moveit.collision_result_right);

    // 返回碰撞检测结果
    return robot_moveit.collision_result_right.collision;
}

/**
 * 对左臂进行碰撞检测
 * @param p : 左臂各关节位置
 * @return  : 碰撞结果
 */
bool Planner::collisionCheckLeft(VectorXd p)
{
    // 将要检测的角度位置更新到模型中去
    std::vector<double> joint_pos(6, 0);

    for(unsigned i = 0; i < 6; i++)
        joint_pos[i] = p[i];

    check_collision_state.setJointGroupPositions(robot_moveit.joint_model_group_left, joint_pos);

    // 清除原来计算结果
    robot_moveit.collision_result_left.clear();

    // 碰撞检测
    robot_moveit.planning_scene->checkCollision(*robot_moveit.collision_request_left,
                                                 robot_moveit.collision_result_left);

    // 返回碰撞检测结果
    return robot_moveit.collision_result_left.collision;
}

void Planner::publishJointState(VectorXd jntLeft, VectorXd jntRight)
{
/*    // 建立jointstate变量，用来发布关节信息
    JointState joint = InitJointState(JntArray(16));

    // 设置关节位置值
    for(unsigned int i = 0; i < 6; i++)
    {
        joint.position[i] = jntLeft[i];
        joint.position[i + 8] = jntRight[i];
    }

    joint.header.stamp = ros::Time::now();

    pub->publish(joint);
    */
}
