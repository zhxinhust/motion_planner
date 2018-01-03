#ifndef _DUAL_ARM_ROBOT
#define _DUAL_ARM_ROBOT

#include <kdl/frames.hpp>
#include <chain.hpp>
#include <frames_io.hpp>
#include <trajectory.hpp>
#include <trajectory_segment.hpp>
#include <trajectory_stationary.hpp>
#include <trajectory_composite.hpp>
#include <trajectory_composite.hpp>
#include <velocityprofile_trap.hpp>
#include <path_roundedcomposite.hpp>
#include <rotational_interpolation_sa.hpp>
#include <utilities/error.h>
#include <trajectory_composite.hpp>
#include "ros/ros.h"
#include "std_msgs/String.h"
#include <sstream>
#include <sensor_msgs/JointState.h>
#include <math.h>

#include <kdl/chainfksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolver.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/jacobian.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>

#include <trajectory.hpp>
#include <trajectory_segment.hpp>
#include <trajectory_stationary.hpp>
#include <trajectory_composite.hpp>
#include <trajectory_composite.hpp>
#include <velocityprofile_trap.hpp>

#include <path_roundedcomposite.hpp>
#include <rotational_interpolation_sa.hpp>
#include <utilities/error.h>
#include <trajectory_composite.hpp>

#include <path_line.hpp>
#include <path_circle.hpp>

#include <trac_ik/trac_ik.hpp>

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

#include "mymath.hpp"
//#include "fgmm++.hpp"
using namespace KDL;

class Joint_limit
{
public:
    double lower;
    double upper;

    Joint_limit();
    ~Joint_limit();

    void set_value(double lower, double upper);
};

class Dual_arm_joint_limit
{
public:
    Joint_limit joint[6];

    Dual_arm_joint_limit();
    ~Dual_arm_joint_limit();

    void init_values(double lower[6], double upper[6]);
};
/******************************************************
// Robot 运动学类，里面包含正逆向运动学求解相关的函数及成员
******************************************************/

class Robot_Kinematics_Annalytical
{
public:
    Robot_Kinematics_Annalytical(KDL::Frame baseFrame_left, KDL::Frame baseFrame_right, KDL::Frame toolFrame_left, KDL::Frame toolFrame_rigth);
    Robot_Kinematics_Annalytical();

    void FK_left(KDL::JntArray jnt_in, KDL::Frame & frame_out);

    // 右臂的正向运动学
    void FK_right(KDL::JntArray jnt_in, KDL::Frame & frame_out);

    void IK_analytical_left(JntArray jnt_init, Frame frame, JntArray &jnt_out);

    void IK_analytical_right(JntArray jnt_init, Frame frame, JntArray &jnt_out);

protected:
    KDL::Frame baseFrame_left_;  // 从全局坐标系到左臂基座标系的变换矩阵
    KDL::Frame baseFrame_right_; // 从全局坐标系到右臂基座标系的变换矩阵
    KDL::Frame toolFrame_left_;  // 从机器人末端到指定tool坐标系的变换矩阵
    KDL::Frame toolFrame_right_;
    int IKSelection(double *q, int qnum, KDL::JntArray jnt_init, KDL::JntArray &jnt_out);
};


class Robot_KIN
{
public:
    TRAC_IK::TRAC_IK *arm_right;    // 右臂Trac ik求解器
    TRAC_IK::TRAC_IK *arm_left;     // 左臂Trac IK求解器
    unsigned int nrOfJoints_right;           // 右臂关节数量
    unsigned int nrOfJoints_left;            // 左臂关节数量
    Dual_arm_joint_limit joint_limits_left, joint_limits_right;
    JntArray jnt_init_left;     // 求逆解时需要用到的初始值
    JntArray jnt_init_right;
    JntArray jnt_current_left;
    JntArray jnt_current_right;

    //Gmm Gmm_dual;
    //Gmm Gmm_left;
    //Gmm Gmm_right;
private:
    Chain chain_left;       // 左臂运动链
    Chain chain_right;      // 右臂运动链
    ChainFkSolverPos_recursive *fk_solver_left;  // 左臂正向运动学求解器
    ChainFkSolverPos_recursive *fk_solver_right; // 右臂正向运动学求解器

    ChainJntToJacSolver *jac_solver_left;   // 左臂雅克比求解器
    ChainJntToJacSolver *jac_solver_right;  // 右臂雅克比求解器
    ChainIkSolverVel_pinv *ik_vel_left;     // 左臂逆向速度求解器
    ChainIkSolverVel_pinv *iK_vel_right;

    double a2, a3, d1, d4, d5, d6;


    Eigen::Vector4d w0, w1, w2, w3, w4, w5;    // 各轴方向
    Eigen::Vector4d p0, p1, p2, p3, p4, p5;    // 各轴上的点

    Eigen::Matrix4d frame_base_left, frame_base_right; // support到各臂基座标的变换矩阵
    Eigen::Matrix4d gst;   // 在各关节值为0时，末端坐标值

    int getUr3Jac(JntArray q_in, Eigen::Matrix4d frame_base, Jacobian &jac);

    const robot_state::JointModelGroup *joint_model_group_left;
    const robot_state::JointModelGroup *joint_model_group_right;

public:
    // 构造函数
    Robot_KIN(const std::string &_urdf_param, const std::string &_chain_start_left, const std::string &_chain_end_left, const std::string &_chain_start_right, const std::string &_chain_end_right, double _time_out = 0.005);

    //
    ~Robot_KIN();

    // 左臂的正向运动学
    void FK_left(JntArray jnt_in, Frame & frame_out);

    // 右臂的正向运动学
    void FK_right(JntArray jnt_in, Frame & frame_out);

    // 左臂的逆向运动学
    int IK_left(Frame frame_in, JntArray &jnt_out, JntArray init_jnt = JntArray(6));

    // 右臂的逆向运动学
    int IK_right(Frame frame_in, JntArray &jnt_out, JntArray init_jnt = JntArray(6));

    int JntToJac_left(const JntArray& q_in, Jacobian& jac, int segmentNR=-1);

    int JntToJac_right(const JntArray& q_in, Jacobian& jac, int segmentNR=-1);

    int IK_vel_left (const JntArray &q_in, const Twist &v_in, JntArray &qdot_out);

    int IK_vel_right (const JntArray &q_in, const Twist &v_in, JntArray &qdot_out);

    int IK_analytical_left(Frame frame, double *q_sol, double q6);

    int IK_analytical_right(Frame frame, double *q_sol, double q6);

    int IK_analytical_left(JntArray jnt_init, Frame frame, JntArray &jnt_out);

    void IK_analytical_left(Frame frame, JntArray &jnt_out);

    int IK_analytical_right(JntArray jnt_init, Frame frame, JntArray &jnt_out);

    void IK_analytical_right(Frame frame, JntArray &jnt_out);

    int GetJac_left(JntArray q_in, Jacobian &jac);

    int GetJac_right(JntArray q_in, Jacobian &jac);

    Eigen::VectorXd Get_singularvalues_left(JntArray q_in);

    Eigen::VectorXd Get_singularvalues_right(JntArray q_in);

    void setInitJnt(JntArray jnt_left, JntArray jnt_right);

    void setInitJntAsCurrent();

    void setCurrentJnt(JntArray jnt_left, JntArray jnt_right);
};

/**
 * 将Frame中的数据拷贝到pose中
 * @param frame_in
 * @param pose_out
 */
void frame2pose(Frame frame_in, geometry_msgs::Pose &pose_out);

/**
 * 数据转换，将pose中的数据拷贝到Frame中去
 * @param pose_in      输入的pose
 * @param frame_out    输出的frame
 */
void pose2frame(geometry_msgs::Pose pose_in, Frame &frame_out);

/**
 * 数据类型转换，将JntArray中的数据拷贝到std::vector<double>中去
 * @param jntarray_in      输入的JntArray
 * @param jntvector_out    输出的vector
 */
void jntArray2jntVector(JntArray jntarray_in, std::vector<double> &jntvector_out);

/**
* 将std::vector<double>中的数据拷贝到JntArray中去
* @param jntvector_in     输入的vector
* @param jntArray_out     输出的JntArray
*/
void jntVector2jntArray(std::vector<double> jntvector_in, JntArray &jntArray_out);

/**
* 初始化JointState变量，用来发布JointState消息，驱动Rviz机器人运动
* @param initjnt ： 初始化值
* @return ： 生成的JointState变量
*/
sensor_msgs::JointState InitJointState(KDL::JntArray initjnt);

/**
 * 设置JointState中的关节值
 * @param jnt_left  ：左臂关节值
 * @param jnt_right ：右臂关节值
 * @param joint     ：JointState变量
 */
void setJointState(JntArray jnt_left, JntArray jnt_right, sensor_msgs::JointState &joint);

void setJointState(JntArray jnt_pos_left, JntArray jnt_pos_right, JntArray jnt_vel_left, JntArray jnt_vel_right, sensor_msgs::JointState &joint);

void copyDataFromArmFrame(Frame frame_arm, float *data);

void getArmFrameFromData(float *data, Frame & frame_left, Frame & frame_right);
void getArmFrameFromVector(Eigen::VectorXd vector_in, Frame &frame_left, Frame &frame_right);
void getArmFrameFromCenter(Frame frame_center, float r, Frame &frame_left, Frame &frame_right);
void getCenterFrameFromData(float *data, Frame & frame_center);
void getVelFromCenter(Frame frame_center, float r, Twist vel_center, Twist &vel_left, Twist &vel_right);

void random2vectorxd(float *random_pos_in, Eigen::VectorXd & vector_out, int dim);
void vectorxd2random(Eigen::VectorXd vector_in, float *random_pos_out, int dim);
void getCenterFrameFromVector(Eigen::VectorXd vector, Frame & frame_center);

#endif
