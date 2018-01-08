//
// Created by zhxin on 17-8-4.
//

#include "RRTPlanner.h"

//#include "RRT_path_planning.h"
#include <fstream>

#define PI 3.1415926

int findAccessFileName();

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

    JntArray str_jnt(6), goal_jnt(6);
    str_jnt.data << -1.69642, -1.51015,  2.20146, 2.45029, 0.125622, -2.35619;
    goal_jnt.data <<-0.448798, -0.688315,   1.43194,   2.39796,    -1.122,  -2.35619;
    
    Robot_moveit robotMoveit;
    robot_state::RobotState &check_collision_state = robotMoveit.planning_scene->getCurrentStateNonConst();

    // 建立RRT规划器类
    RRTPlanner rrtPlanner(check_collision_state, robotMoveit);

    // 初始化RRT规划器
    rrtPlanner.Init(str_jnt.data, goal_jnt.data);

    // 进行规划
    rrtPlanner.plan();

    // 显示规划运动
    VectorXd str, goal, v_direct, v_dis_temp, v_temp;
    double ds = 0.05, s = 0;

    int index = findAccessFileName();

    char *filename = new char[100];
    char num[5]={0};

    // 保存文件
    strcpy(filename, "/home/zhaoxin/Code/catkin_gmm_multirrt/src/motion_planners/data/samples");
    sprintf(num, "%d", index);
    strcat(filename, num);
    strcat(filename,".txt");
    std::ofstream samples_file(filename);

    strcpy(filename, "/home/zhaoxin/Code/catkin_gmm_multirrt/src/motion_planners/data/path");
    strcat(filename, num);
    strcat(filename,".txt");
    std::ofstream path_file(filename);

    for(int i = 0; i < rrtPlanner.tree_start.size(); i++)
        samples_file << rrtPlanner.tree_start[i].vector.transpose() << '\n';

    for(int i = 0; i < rrtPlanner.tree_goal.size(); i++)
        samples_file << rrtPlanner.tree_goal[i].vector.transpose() << '\n';

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
            {
                joint.position[i + 9] = v_temp[i];
            }
            path_file << v_temp.transpose() << '\n';
//            outfile << '\n';

            // 发布关节位置
            joint.header.stamp = ros::Time::now();
            pub.publish(joint);


            loop_rate.sleep();
        }
    }
}

using  namespace std;
int findAccessFileName()
{
    FILE *_file;
    std::fstream fs;

    char *filename;
    filename = new char[100];

    char num[5]={0};
    int i = 0;
    do{
        strcpy(filename, "/home/zhaoxin/Code/catkin_gmm_multirrt/src/motion_planners/data/samples");
        sprintf(num, "%d", i);
        strcat(filename, num);
        strcat(filename,".txt");
        _file = fopen(filename,"r");
//        fs.open(filename, ios::in);
        i++;
    }while(_file != NULL);

    return i - 1;

}