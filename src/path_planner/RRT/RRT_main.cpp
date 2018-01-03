//
// Created by zhxin on 17-8-4.
//

#include "RRTPlanner.h"

//#include "RRT_path_planning.h"

#define PI 3.1415926

int main(int argc, char **argv) {
    ros::init(argc, argv, "dual_arm_motion_planning");
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
    Robot_KIN rob_kin(urdf_param, "world", "left_tool0", "world", "right_tool0", 0.005);

    // 构建moveit所需要的类
   // Robot_moveit *robot_moveit = new Robot_moveit;

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

//    for(int i = 0; i < 6; i++)
//        joint.position[i + 9] = jnt_str_pos_right.data[i];
//    joint.header.stamp = ros::Time::now();
//    pub.publish(joint);
//
//    for(int i = 0; i < 6; i++)
//        joint.position[i + 9] = jnt_ik.data[i];
//    joint.header.stamp = ros::Time::now();
//    pub.publish(joint);

    Robot_moveit robotMoveit;
    robot_state::RobotState &check_collision_state = robotMoveit.planning_scene->getCurrentStateNonConst();

    // 建立RRT规划器类
    RRTPlanner rrtPlanner(check_collision_state, robotMoveit);

    // 初始化RRT规划器
    rrtPlanner.Init(jnt_str_pos_right.data, jnt_ik.data);

    // 进行规划
    rrtPlanner.plan();

    // 显示规划运动
    VectorXd str, goal, v_direct, v_dis_temp, v_temp;
    double ds = 0.05, s = 0;

    for(int i = 0; i < (int)rrtPlanner.path.size() - 1; i++)
    {
        str = rrtPlanner.path[i];
        goal = rrtPlanner.path[i + 1];
        v_dis_temp = goal - str;
        v_direct = v_dis_temp / v_dis_temp.norm();
        s = 0;
        while(s < v_dis_temp.norm())
        {
            v_temp = str + s * v_direct;
            s = s + ds;

            for(unsigned int i = 0; i < 6; i++)
                joint.position[i + 9] = v_temp[i];

            // 发布关节位置
            joint.header.stamp = ros::Time::now();
            pub.publish(joint);

            loop_rate.sleep();
        }
    }
}