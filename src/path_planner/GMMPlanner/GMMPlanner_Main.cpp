//
// Created by zhxin on 17-8-24.
//

#include "GMMGuidedPlanner.hpp"

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

    sleep(3.0);

    ROS_INFO("waiting for Rviz, please wait...");

    JointState joint = InitJointState(JntArray(17));

    std::vector<VectorXd> rrt_path;

    VectorXd p_str, p_goal;

    p_str.resize(6);
    p_goal.resize(6);

    srand((unsigned)time(0));   // 生成种子

//    Frame strFrame = Frame(Rotation(1, 0, 0, 0, 0, -1, 0, 1, 0), Vector(-0.4, 0.1, 0.2));
    Frame strFrame = Frame(KDL::Rotation::Quaternion(-0.5, 0.5, -0.5, 0.5), Vector(0.15, -0.24, 1.046));
    Frame goalFrame = Frame(KDL::Rotation::Quaternion(-0.5, 0.5, -0.5, 0.5), Vector(0.686, -0.24, 1.286));

    JntArray strJntRight(6), goalJntRight(6);

    rob_kin.IK_right(strFrame, strJntRight, JntArray(6));
    rob_kin.IK_right(goalFrame, goalJntRight, strJntRight);

    std::cout << "start joint:" << strJntRight.data.transpose() << std::endl;
    std::cout << "goal joint" << goalJntRight.data.transpose() << std::endl;

    for(unsigned int i = 0; i < 6; i++)
    {
//        joint.position[i] = 0;//0～5左臂，6~8左爪
        joint.position[i + 9] = strJntRight.data[i];//9~14右臂，15～16右爪
    }
    // 发布关节位置
    joint.header.stamp = ros::Time::now();
    pub.publish(joint);

    for(unsigned int i = 0; i < 6; i++)
    {
//        joint.position[i] = 0;//0～5左臂，6~8左爪
        joint.position[i + 9] = goalJntRight.data[i];//9~14右臂，15～16右爪
    }
    // 发布关节位置
    joint.header.stamp = ros::Time::now();
    pub.publish(joint);

    Robot_moveit robotMoveit;
    robot_state::RobotState &check_collision_state = robotMoveit.planning_scene->getCurrentStateNonConst();

    // 建立规划器类
    GMMGuidedPlanner gmmplanner(check_collision_state, robotMoveit);

    // 初始化规划器
    gmmplanner.Init(strJntRight.data, goalJntRight.data);

    // 进行规划
    gmmplanner.constructRoadMap();

    gmmplanner.plan();

    // 显示规划运动
    VectorXd str, goal, v_direct, v_dis_temp, v_temp;
    double ds = 0.015, s = 0;

//    FILE *wfile = fopen("/home/zhxin/catkin_ws/src/dual_arm_robot/src/testdata/data.txt", "w");

    for(int i = 0; i < (int)gmmplanner.path.size() - 1; i++)
    {
        str = gmmplanner.path[i];
        goal = gmmplanner.path[i + 1];
        v_dis_temp = goal - str;
        v_direct = v_dis_temp / v_dis_temp.norm();
        s = 0;
        while(s < v_dis_temp.norm())
        {
            v_temp = str + s * v_direct;
            s = s + ds;

//            std::cout << v_temp << std::endl;
            for(unsigned int i = 0; i < 6; i++)
            {
                joint.position[i + 9] = v_temp[i];
//                fprintf(wfile, "%lf, ", v_temp[i]);
            }
//            fprintf(wfile, "\n");


            // 发布关节位置
            joint.header.stamp = ros::Time::now();
            pub.publish(joint);

            loop_rate.sleep();
        }
    }
//    fclose(wfile);

    return 1;
}
