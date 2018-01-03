
#include "RRTStarPlanner.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <visualization_msgs/Marker.h>
#include <cmath>
#define PI 3.1415926

char *findAccessFileName();

int main(int argc, char **argv) {
    ros::init(argc, argv, "dual_arm_assembly_bottle");
    ros::NodeHandle node_handle;
    ros::AsyncSpinner spinner(1);
    spinner.start();

    // 建立发布器
    ros::Publisher pub = node_handle.advertise<sensor_msgs::JointState>("/joint_states", 10);

    // 通信频率
    ros::Rate loop_rate(20);

    std::string urdf_param;

    double timeout;
    node_handle.param("timeout", timeout, 0.005);
    node_handle.param("urdf_param", urdf_param, std::string("/robot_description"));

    // 建立正逆向运动学模型
    Robot_KIN rob_kin(urdf_param, "uu_support", "left_tool0", "uu_support", "right_tool0", 0.005);

    // 构建moveit所需要的类
    // Robot_moveit *robot_moveit = new Robot_moveit;

    sleep(5.0);

    ROS_INFO("waiting for Rviz, please wait...");

    JointState joint = InitJointState(JntArray(17));

    //使用规划器
    std::vector<VectorXd> rrt_path;

    VectorXd p_str, p_goal;

    p_str.resize(6);
    p_goal.resize(6);

    Frame strFrame = Frame(Rotation(1, 0, 0, 0, 0, -1, 0, 1, 0), Vector(-0.4, 0.1, 0.2));
    Frame goalFrame = Frame(Rotation(1, 0, 0, 0, 0, -1, 0, 1, 0), Vector(-0.15, -0.14, 0.4));

    JntArray strJnt(6);
    JntArray goalJnt(6);

    rob_kin.IK_right(strFrame, strJnt, JntArray(6));
    rob_kin.IK_right(goalFrame, goalJnt, strJnt);

    p_str = strJnt.data;
    p_goal = goalJnt.data;

    std::cout << "start Joint is: " << p_str.transpose() << std::endl;
    std::cout << "end joint is: " << p_goal.transpose() << std::endl;

    char *filename = findAccessFileName();
    std::cout << filename << std::endl;



    for(unsigned int i = 0; i < 6; i++)
    {
        joint.position[i] = 0;//0～5左臂，6~8左爪
        joint.position[i + 9] = p_str[i];//9~14右臂，15～16右爪
    }
    // 发布关节位置
    joint.header.stamp = ros::Time::now();
    pub.publish(joint);

    for(unsigned int i = 0; i < 6; i++)
    {
        joint.position[i + 9] = p_goal[i];//9~14右臂，15～16右爪
    }
    // 发布关节位置
    joint.header.stamp = ros::Time::now();
    pub.publish(joint);

    Robot_moveit robotMoveit;
    robot_state::RobotState &check_collision_state = robotMoveit.planning_scene->getCurrentStateNonConst();

    // 建立RRT规划器类
    //RRTStarPlanner rrtStarPlanner(check_collision_state, robotMoveit);
    RRTPlanner rrtStarPlanner(check_collision_state, robotMoveit);

    // 初始化RRT规划器
    rrtStarPlanner.Init(p_str, p_goal);

    // 进行规划
    rrtStarPlanner.plan();

    std::vector<VectorXd> optimalPath_ = rrtStarPlanner.path;
    int row = (int)rrtStarPlanner.path.size();
//    std::vector<VectorXd> optimalPath_;
//    int row = motionPlanning(optimalPath_);


    // 显示规划运动
    VectorXd str, goal, v_direct, v_dis_temp, v_temp;
    double ds = 0.08, s = 0;


    FILE *wfile = fopen(filename, "w");

    for(int i = 0; i < row - 1; i++)
    {

        //使用规划器
        str = optimalPath_[i];
        goal = optimalPath_[i + 1];

        v_dis_temp = goal - str;
        v_direct = v_dis_temp / v_dis_temp.norm();
        s = 0;
        while(s < v_dis_temp.norm())
        {
            v_temp = str + s * v_direct;
            s = s + ds;

            for(unsigned int i = 0; i < 6; i++)
            {
                //joint.position[i] = v_temp[i];//0～5左臂，6~8左爪
                joint.position[i + 9] = v_temp[i];//9~14右臂，15～16右爪
                fprintf(wfile, "%f ", v_temp[i]);
            }
            fprintf(wfile, "\n");

            // 发布关节位置
            joint.header.stamp = ros::Time::now();
            pub.publish(joint);

            loop_rate.sleep();
        }
    }
    fclose(wfile);
}


using  namespace std;
char *findAccessFileName()
{
    FILE *_file;

    char *filename;
    filename = new char[30];

    char num[5]={0};
    int i = 0;

    do{
        strcpy(filename, "/home/zhaoxin/catkin_ws/src/dual_arm_robot/src/path_planning/data/data");
        sprintf(num, "%d", i);
        strcat(filename, num);
        strcat(filename,".txt");
        _file = fopen(filename,"r");
        i++;
    }while(_file != NULL);

    return filename;

}