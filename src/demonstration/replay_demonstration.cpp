//
// Created by zhaoxin on 18-1-4.
//

#include "dual_arm_robot.hpp"
#include <fstream>

using namespace std;
int main(int argc, char **argv) {
    ros::init(argc, argv, "dual_arm_motion_planning");
    ros::NodeHandle node_handle;
    ros::AsyncSpinner spinner(1);
    spinner.start();

    // 建立发布器
    ros::Publisher pub = node_handle.advertise<sensor_msgs::JointState>("/joint_states", 150);

    // 通信频率
    ros::Rate loop_rate(125);

    sensor_msgs::JointState joint = InitJointState(JntArray(17));

    ifstream infile;

    infile.open("/home/zhaoxin/Code/catkin_gmm_multirrt/src/motion_planners/data/demonstration/data4.txt");

    assert(infile.is_open());   //若失败,则输出错误消息,并终止程序运行

    string s;
    char buffer[256] = {0};
    getline(infile,s);

    double a;
    while (! infile.eof() )
    {
        infile.getline (buffer,150);
//        std::cout << buffer <<std::endl;
        sscanf(buffer,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&a,&joint.position[9],&joint.position[10], &joint.position[11], &joint.position[12], &joint.position[13], &joint.position[14]);
        joint.header.stamp = ros::Time::now();
        pub.publish(joint);
        loop_rate.sleep();
    }
}