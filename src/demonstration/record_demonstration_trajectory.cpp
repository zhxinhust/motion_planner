//
// Created by zhaoxin on 18-1-4.
//

#include <sensor_msgs/JointState.h>
#include <fstream>
#include "ros/ros.h"

using  namespace std;
using  namespace sensor_msgs;

int findAccessFileName()
{
    FILE *_file;
    std::fstream fs;

    char *filename;
    filename = new char[100];

    char num[5]={0};
    int i = 0;
    do{
        strcpy(filename, "/home/zhaoxin/Code/catkin_gmm_multirrt/src/motion_planners/data/path");
        sprintf(num, "%d", i);
        strcat(filename, num);
        strcat(filename,".txt");
        _file = fopen(filename,"r");
//        fs.open(filename, ios::in);
        i++;
    }while(_file != NULL);

    return i - 1;

}


std::ofstream samples_file;

FILE *wf;

int datanum = 0;
void record_file(const sensor_msgs::JointState &msg)
{
    datanum++;
    if(datanum % 4 == 0){
        std::cout << datanum << std::endl;
        for(int i = 0; i < 6; i++)
            fprintf(wf, "%f ", msg.position[i]);
        samples_file << '\n';
    }

}

int main(int argc, char **argv) {
    ros::init(argc, argv, "dual_arm_motion_planning");
    ros::NodeHandle node_handle;
    ros::AsyncSpinner spinner(1);
    spinner.start();

    int index = findAccessFileName();
    char *filename = new char[100];
    char num[5]={0};

    // 保存文件
    strcpy(filename, "/home/zhaoxin/Code/catkin_gmm_multirrt/src/motion_planners/data/path");
    sprintf(num, "%d", index);
    strcat(filename, num);
    strcat(filename,".txt");
    wf = fopen(filename, "w");
//    samples_file.open(filename);

    ros::Subscriber sub = node_handle.subscribe("joint_states", 1000, record_file);

    ros::spin();

    return 0;
}
