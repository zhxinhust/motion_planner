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

    // 构建moveit所需要的类
    // Robot_moveit *robot_moveit = new Robot_moveit;

    sleep(1.0);


}