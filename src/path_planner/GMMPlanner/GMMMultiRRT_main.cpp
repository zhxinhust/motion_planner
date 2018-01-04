//
// Created by zhaoxin on 17-12-28.
//

#include "GMMGuidedMultiRRTPlanner.hpp"

int main(int argc, char **argv)
{
    ros::init(argc, argv, "dual_arm_assembly_bottle");
    ros::NodeHandle node_handle;
    ros::AsyncSpinner spinner(1);
    spinner.start();

    // 建立发布器
    ros::Publisher pub = node_handle.advertise<sensor_msgs::JointState>("/joint_states", 10);

    // 通信频率
    ros::Rate loop_rate(40);

    std::string urdf_param;

    double timeout;
    node_handle.param("timeout", timeout, 0.005);
    node_handle.param("urdf_param", urdf_param, std::string("/robot_description"));

    // 建立正逆向运动学模型
    Robot_KIN rob_kin(urdf_param, "uu_support", "left_tool0", "uu_support", "right_tool0", 0.005);

    sleep(1.0);

    ROS_INFO("waiting for Rviz, please wait...");

    JointState joint = InitJointState(JntArray(17));

    std::vector<VectorXd> rrt_path;

    VectorXd p_str, p_goal;

    p_str.resize(6);
    p_goal.resize(6);

    JntArray jnt_init_right(6), jnt_str_pos_right(6), jnt_str_pos_left(6), jnt_ik(6);
    jnt_str_pos_left.data << -0.0706, -0.8472, -2.329, -0.8472, -3.8828, 0;
    jnt_str_pos_right.data << 0.0706, -1.9767, 2.1179, -1.2707, -2.6827, 0;
    jnt_init_right.data << -0.38453094036000035, -0.9826901809199997, 1.2805131641400003, -0.25635396023999935, 1.1963184811199996, 0.7263362206799995;

    Frame strFrame = Frame(KDL::Rotation::Quaternion(-0.5, 0.5, -0.5, 0.5), Vector(0.15, -0.24, 1.046));
    Frame goalFrame = Frame(KDL::Rotation::Quaternion(-0.5, 0.5, -0.5, 0.5), KDL::Vector(0.686, -0.24, 1.286));

    rob_kin.IK_right(strFrame, jnt_str_pos_right, JntArray(6));
    rob_kin.IK_right(goalFrame, jnt_ik, jnt_str_pos_right);

    Robot_moveit robotMoveit;
    robot_state::RobotState &check_collision_state = robotMoveit.planning_scene->getCurrentStateNonConst();

    GMMGuidedMultiRRTPlanner gmmMultirrtPlanner(check_collision_state, robotMoveit);

    // 初始化RRT规划器
    gmmMultirrtPlanner.Init(jnt_str_pos_right.data, jnt_ik.data);

    // 进行规划
    gmmMultirrtPlanner.plan();
}