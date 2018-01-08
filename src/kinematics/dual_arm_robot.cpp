/**
 *
 */
#include "dual_arm_robot.hpp"
#include "ur_kin.h"
#include "mymath.hpp"

using namespace KDL;
using namespace sensor_msgs;
using namespace ur_kinematics;
using namespace Eigen;


inline int IK_analytical(Frame frame_base, Frame frame_obj, double *q_sol, double q6);
inline int IK_selection(JntArray jnt_init, double *q_sol, int q_num, JntArray &jnt_out);


inline void copyDataFromFrameToArray(KDL::Frame frame_obj, double *T)
{
    for(int i = 0; i < 3; i++)
    {
        // 获取旋转矩阵
        *(T + i) = frame_obj.M.data[i];
        *(T + 4 + i) = frame_obj.M.data[3 + i];
        *(T + 8 + i) = frame_obj.M.data[6 + i];

        *(T + 4*(i + 1) - 1) = frame_obj.p.data[i];     // 获取位置量
    }
    *(T + 15) = 1;
}

inline void copyDataFromArrayToFrame(double *T, KDL::Frame &frame_obj)
{
    for(int i = 0; i < 3; i++)
    {
        // 获取旋转矩阵
        frame_obj.M.data[i] = *(T + i);
        frame_obj.M.data[3 + i] = *(T + 4 + i);
        frame_obj.M.data[6 + i] = *(T + 8 + i);

        frame_obj.p.data[i] = *(T + 4*(i + 1) - 1);     // 获取位置量
    }
}

inline void copyDataFromJntArrayToArray(KDL::JntArray jnt, double *q)
{
    for(int i = 0; i < 6; i++)
    {
        *(q + i) = jnt.data(i);
    }
}

inline void copyDataFromArrayToJntArray(double *q, KDL::JntArray &jnt_out)
{
    for(int i = 0; i < 6; i++)
    {
        jnt_out.data[i] = *(q + i);
    }
}

/**
 * 设置关节运动范围
 * @param lower_
 * @param upper_
 */
void Joint_limit::set_value(double lower_, double upper_)
{
    lower = lower_;
    upper = upper_;
}

Joint_limit::Joint_limit() {}

Joint_limit::~Joint_limit() {}

Dual_arm_joint_limit::Dual_arm_joint_limit() {}

Dual_arm_joint_limit::~Dual_arm_joint_limit() {};

/**
 * 初始化单个机械臂
 * @param lower_
 * @param upper_
 */
void Dual_arm_joint_limit::init_values(double *lower_, double *upper_)
{

    for(int i = 0; i < 6; i++)
    {
        joint[i].upper = upper_[i];
        joint[i].lower = lower_[i];
    }
}


/**
 * 双臂机器人正逆向运动学类构造函数
 * @param _urdf_param           URDF文件名
 * @param _chain_start_left     左臂运动链起始连杆名称
 * @param _chain_end_left       左臂运动链终止连杆名称
 * @param _chain_start_right    右臂运动链起始连杆名称
 * @param _chain_end_right      右臂运动链终止连杆名称
 * @param _time_out             逆解求解器超时时间
 */
Robot_KIN::Robot_KIN(const std::string &_urdf_param,
                     const std::string &_chain_start_left, const std::string &_chain_end_left,
                     const std::string &_chain_start_right, const std::string &_chain_end_right,
                     double _time_out)
{
    double eps = 1e-5;  // 求解精度
    // 构造右臂的TRAC IK求解器
    arm_left = new TRAC_IK::TRAC_IK(_chain_start_left, _chain_end_left, _urdf_param, _time_out, eps);

    // 构造左臂的TRAC IK求解器
    arm_right = new TRAC_IK::TRAC_IK(_chain_start_right, _chain_end_right, _urdf_param, _time_out, eps);

    ROS_INFO("Trac IK finished!");

    arm_left->getKDLChain(chain_left);      // 获取左臂运动学链
    arm_right->getKDLChain(chain_right);    // 获取右臂运动学链

    nrOfJoints_left = chain_left.getNrOfJoints();   // 保存左臂关节数
    nrOfJoints_right = chain_right.getNrOfJoints(); // 保存右臂关节数

    jac_solver_left = new ChainJntToJacSolver(chain_left);
    jac_solver_right = new ChainJntToJacSolver(chain_right);

    fk_solver_left = new ChainFkSolverPos_recursive(chain_left);    // 建立左臂正向运动学求解器
    fk_solver_right = new ChainFkSolverPos_recursive(chain_right);    // 建立右臂正向运动学求解器

    ik_vel_left = new ChainIkSolverVel_pinv(chain_left);     // 左臂逆向速度求解器
    iK_vel_right = new ChainIkSolverVel_pinv(chain_right);   // 右臂逆向速度求解器

    jnt_init_left = JntArray(nrOfJoints_left);
    jnt_init_right = JntArray(nrOfJoints_right);

    jnt_current_left = JntArray(nrOfJoints_left);
    jnt_current_right = JntArray(nrOfJoints_right);

    // 连杆参数值
    a2 = 0.24365;
    a3 = 0.21325;

    d1 = 0.1519;
    d4 = 0.11235;
    d5 = 0.08535;
    d6 = 0.0819;

    w0 << 0.0, 0.0, 1.0, 0.0;
    w1 << 0, -1, 0, 0;
    w2 << 0, -1, 0, 0;
    w3 << 0, -1, 0, 0;
    w4 << 0, 0, -1, 0;
    w5 << 0, -1, 0, 0;

    p0 << 0, 0, 0, 1;
    p1 << 0, 0, d1, 1;
    p2 << -a2, 0, d1, 1;
    p3 << -a2 - a3, 0, d1, 1;
    p4 << -a2 - a3, -d4, d1, 1;
    p5 << -a2 - a3, -d4, d1 - d5, 1;

    gst << 1, 0, 0, -a2 - a3,
         0, 0, -1, -d4 - d6,
         0, 1, 0, d1 - d5,
         0, 0, 0, 1;

    // 由世界坐标系到左臂基座标系的变换矩阵
    Frame world_frame_left = Frame(Rotation::RPY(2.526113, -0.523599, 2.526113), Vector(-0.001, 0.195, -0.015));

    Rotation M = world_frame_left.M;

    frame_base_left << M.data[0], M.data[1], M.data[2], -0.001,
                       M.data[3], M.data[4], M.data[5], 0.195,
                       M.data[6], M.data[7], M.data[8], -0.015,
                       0, 0, 0, 1;

    // 由世界坐标系到右臂基座标系的变换矩阵
    Frame world_frame_right = Frame(Rotation::RPY(2.526113, 0.523599, 0.615480), Vector(-0.001, -0.195, -0.015));

    M = world_frame_right.M;

    frame_base_right<<M.data[0], M.data[1], M.data[2], -0.001,
                      M.data[3], M.data[4], M.data[5], -0.195,
                      M.data[6], M.data[7], M.data[8], -0.015,
                      0, 0, 0, 1;

    // 设定各关节运动极限
    double lower_value[6] = {-2.0 * M_PI};
    double upper_value[6] = {2.0 * M_PI};

    joint_limits_left.init_values(lower_value, upper_value);
    joint_limits_right.init_values(lower_value, upper_value);

}

/**
 * 析构函数
 */
Robot_KIN::~Robot_KIN()
{
    free(jac_solver_left);
    free(jac_solver_right);
    free(fk_solver_left);
    free(fk_solver_right);
    free(ik_vel_left);
    free(iK_vel_right);
}

/**
 * 左臂的正向运动学
 * @param jnt   :输入的关节值
 * @param frame :输出的笛卡尔空间frame
 */
void Robot_KIN::FK_left(JntArray jnt, Frame & frame)
{
    // ChainFkSolverPos_recursive fk_solver_left = ChainFkSolverPos_recursive(chain_left);    // 建立左臂正向运动学求解器
    fk_solver_left->JntToCart(jnt, frame);
}

/**
 * 右臂的正向运动学
 * @param jnt   输入的关节值
 * @param frame 输出的笛卡尔空间frame
 */
void Robot_KIN::FK_right(JntArray jnt, Frame & frame)
{
    //ChainFkSolverPos_recursive fk_solver_right = ChainFkSolverPos_recursive(chain_right);    // 建立右臂正向运动学求解器
    fk_solver_right->JntToCart(jnt, frame);
}

/**
 * 左臂的逆向运动学
 * @param frame_in  :输入的笛卡尔坐标系下的Frame
 * @param jnt_out   :输出的关节位置
 * @param jnt_init  :初始关节位置，此值默认参数为0
 * @return :如果计算得到逆解则返回1，否则返回其他值
 */
int Robot_KIN::IK_left(Frame frame_in, JntArray &jnt_out, JntArray jnt_init)
{
    return arm_left->CartToJnt(jnt_init, frame_in, jnt_out);
}

/**
 * 右臂的逆向运动学
 * @param frame     笛卡尔坐标
 * @param jnt_out   输出的关节值
 * @param jnt_init  数值求解的初始值，此值默认全为0
 * @return 如果计算得到逆解则返回1，否则返回其他值
 */
int Robot_KIN::IK_right(Frame frame_in, JntArray &jnt_out, JntArray jnt_init)
{
    return arm_right->CartToJnt(jnt_init, frame_in, jnt_out);
}

/**
 * 求解雅克比矩阵
 * @param q_in      ： 关节位置
 * @param jac       ： 输出的雅克比矩阵
 * @param segmentNR
 * @return
 */
int Robot_KIN::JntToJac_left(const JntArray& q_in, Jacobian& jac, int segmentNR)
{
    return jac_solver_left->JntToJac(q_in, jac, segmentNR);
}

/**
 * 求解雅克比矩阵
 * @param q_in      ： 关节位置
 * @param jac       ： 输出的雅克比矩阵
 * @param segmentNR
 * @return
 */
int Robot_KIN::JntToJac_right(const JntArray& q_in, Jacobian& jac, int segmentNR)
{
    return jac_solver_right->JntToJac(q_in, jac, segmentNR);
}

/**
 * 求左臂关节速度
 * @param q_in      ： 关节位置
 * @param v_in      ： 末端速度
 * @param qdot_out  ： 关节速度
 * @return
 */
int Robot_KIN::IK_vel_left (const JntArray &q_in, const Twist &v_in, JntArray &qdot_out)
{
    return ik_vel_left->CartToJnt(q_in, v_in, qdot_out);
}

/**
 * 求右臂关节速度
 * @param q_in      ： 关节位置
 * @param v_in      ： 末端速度
 * @param qdot_out  ： 关节速度
 * @return
 */
int Robot_KIN::IK_vel_right (const JntArray &q_in, const Twist &v_in, JntArray &qdot_out)
{
    return iK_vel_right->CartToJnt(q_in, v_in, qdot_out);
}

/**
 * 用ur_kinematics求解析IK
 * @param frame     ： 末端frame
 * @param q_sol     ： 求得的逆解
 * @param q6        :  6关节值，如果求出来奇异，直接使用此值
 * @return          ： 返回解的数量
 */
int Robot_KIN::IK_analytical_left(Frame frame, double *q_sol, double q6)
{
    // 世界坐标系到机器人基座之间的变换
    Frame world_frame_left = Frame(Rotation::RPY(2.526113, -0.523599, 2.526113), Vector(-0.001, 0.195, -0.015));
    return IK_analytical(world_frame_left, frame, q_sol, q6);
}

/**
 * 用ur_kinematics求解析IK
 * @param frame     ： 末端frame
 * @param q_sol     ： 求得的逆解
 * @param q6        :  6关节值，如果求出来奇异，直接使用此值
 * @return          ： 返回解的数量
 */
int Robot_KIN::IK_analytical_right(Frame frame, double *q_sol, double q6)
{
    // 世界坐标系到机器人基座之间的变换
    Frame world_frame_right = Frame(Rotation::RPY(2.526113, 0.523599, 0.615480), Vector(-0.001, -0.195, -0.015));
    return IK_analytical(world_frame_right, frame, q_sol, q6);
}

/**
 * 利用解析IK求解器求解，同时还需要进行选解
 * @param jnt_init  ： 初始关节角，利用此角度进行选解
 * @param frame     ： 末端frame，描述tool0的坐标
 * @param jnt_out   ： 输出关节角
 * @return          ： 如果成功，返回1，无解，返回-1，存在多个解，则返回解的数量
 */
int Robot_KIN::IK_analytical_left(JntArray jnt_init, Frame frame, JntArray &jnt_out)
{
    double q_sol[8 * 6] = {0};    // 关节逆解
    int q_num;              // 解的数量
    q_num = IK_analytical_left(frame, q_sol, jnt_init.data[5]);

    return IK_selection(jnt_init, q_sol, q_num, jnt_out);
}

void Robot_KIN::IK_analytical_left(Frame frame, JntArray &jnt_out)
{
    IK_analytical_left(jnt_init_left, frame, jnt_out);
}

int Robot_KIN::IK_analytical_right(JntArray jnt_init, Frame frame, JntArray &jnt_out)
{
    double q_sol[8 * 6] = {0};    // 关节逆解
    int q_num;              // 解的数量
    q_num = IK_analytical_right(frame, q_sol, jnt_init.data[5]);

    std::cout<<"IK" << std::endl;
    for(int i = 0; i < q_num; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            std::cout<< q_sol[i * 6 + j] << ' ';
        }
        std::cout<<std::endl;
    }

    return IK_selection(jnt_init, q_sol, q_num, jnt_out);
}

void Robot_KIN::IK_analytical_right(Frame frame, JntArray &jnt_out)
{
    IK_analytical_right(jnt_init_right, frame, jnt_out);
}

// 计算单臂的雅克比矩阵
int Robot_KIN::getUr3Jac(JntArray q_in, Eigen::Matrix4d frame_base, Jacobian &jac)
{
    Eigen::Matrix4d A0, A1, A2, A3, A4, A5;

    Eigen::Matrix4d temp;

    Eigen::MatrixXd P(6, 6);

    P.setZero();
    P.block(0, 0, 3, 3) = frame_base.block(0, 0, 3, 3);
    P.block(3, 3, 3, 3) = frame_base.block(0, 0, 3, 3);

    // 求各变换矩阵
    A0 = POE(w0, p0, q_in.data(0));
    A1 = POE(w1, p1, q_in.data(1));
    A2 = POE(w2, p2, q_in.data(2));
    A3 = POE(w3, p3, q_in.data(3));
    A4 = POE(w4, p4, q_in.data(4));
    A5 = POE(w5, p5, q_in.data(5));

    Vector4d p0b, p1b, p2b, p3b, p4b, p5b, pe;  // 各轴上一点坐标
    Vector4d w0b, w1b, w2b, w3b, w4b, w5b;      // 各轴旋转方向

    // 由于建立的x0-y0-z0与模型相比，x0, y0方向相反，这里乘以一个矩阵变换过来
    temp << -1, 0, 0, 0,
            0, -1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    p0b = temp * p0;
    w0b = temp * w0;

    temp = temp * A0;
    p1b = temp * p1;
    w1b = temp * w1;

    temp = temp * A1;
    p2b = temp * p2;
    w2b = temp * w2;

    temp = temp * A2;
    p3b = temp * p3;
    w3b = temp * w3;

    temp = temp * A3;
    p4b = temp * p4;
    w4b = temp * w4;

    temp = temp * A4;
    p5b = temp * p5;
    w5b = temp * w5;

    // 末端点位置
    Vector4d p00;
    p00 << 0, 0, 0, 1;
    pe = temp * A5 * gst * p00;

    // 计算雅克比矩阵各列 [wi x (pe - pi), wi]

    jac.data.block(0, 0, 3, 1) = w0b.head<3>().cross(pe.head<3>() - p0b.head<3>());
    jac.data.block(3, 0, 3, 1) = w0b.head<3>();

    jac.data.block(0, 1, 3, 1) = w1b.head<3>().cross(pe.head<3>() - p1b.head<3>());
    jac.data.block(3, 1, 3, 1) = w1b.head<3>();

    jac.data.block(0, 2, 3, 1) = w2b.head<3>().cross(pe.head<3>() - p2b.head<3>());
    jac.data.block(3, 2, 3, 1) = w2b.head<3>();

    jac.data.block(0, 3, 3, 1) = w3b.head<3>().cross(pe.head<3>() - p3b.head<3>());
    jac.data.block(3, 3, 3, 1) = w3b.head<3>();

    jac.data.block(0, 4, 3, 1) = w4b.head<3>().cross(pe.head<3>() - p4b.head<3>());
    jac.data.block(3, 4, 3, 1) = w4b.head<3>();

    jac.data.block(0, 5, 3, 1) = w5b.head<3>().cross(pe.head<3>() - p5b.head<3>());
    jac.data.block(3, 5, 3, 1) = w5b.head<3>();

    jac.data = jac.data * P;

    return 1;
}

/**
 * 获取左臂的雅克比矩阵
 * @param q_in ：输入的关节值
 * @param jac  ：输出的雅克比矩阵
 * @return
 */
int Robot_KIN::GetJac_left(JntArray q_in, Jacobian &jac)
{
    return getUr3Jac(q_in, frame_base_left, jac);
}

/**
 * 计算右臂的雅克比矩阵
 * @param q_in ：输入的关节值
 * @param jac  ：输出的雅克比矩阵
 * @return
 */
int Robot_KIN::GetJac_right(JntArray q_in, Jacobian &jac)
{
    return getUr3Jac(q_in, frame_base_right, jac);
}

/**
 * 获取左臂的雅克比矩阵的奇异值，用来计算可操作性等指标
 * @param q_in ： 输入关节值
 * @return     ： 奇异值向量
 */
VectorXd Robot_KIN::Get_singularvalues_left(JntArray q_in)
{
    Jacobian jac(6);
    GetJac_left(q_in, jac);

 //   std::cout << jac.data << std::endl;
    // 奇异值分解
    JacobiSVD<MatrixXd> svd(jac.data, ComputeThinU | ComputeThinV);

    //std::cout << svd.singularValues() << std::endl << std::endl;
    VectorXd sv = svd.singularValues();
    return sv;
}

/**
 * 获取右臂的雅克比矩阵的奇异值，用来计算可操作性等指标
 * @param q_in ： 关节值
 * @return     ： 奇异值向量
 */
VectorXd Robot_KIN::Get_singularvalues_right(JntArray q_in)
{
    Jacobian jac(6);
    GetJac_right(q_in, jac);

    // 奇异值分解
    JacobiSVD<MatrixXd> svd(jac.data, ComputeThinU | ComputeThinV);

    return svd.singularValues();
}

/**
 * 设置求逆解需要用到的关节初始值
 * @param jnt_left
 * @param jnt_right
 */
void Robot_KIN::setInitJnt(JntArray jnt_left, JntArray jnt_right)
{
    jnt_init_left = jnt_left;
    jnt_init_right = jnt_right;
}

/**
 * 将初始角度设为当前角度值
 */
void Robot_KIN::setInitJntAsCurrent()
{
    jnt_init_left = jnt_current_left;
    jnt_init_right = jnt_current_right;
}

/**
 * 设置当前角度值
 * @param jnt_left
 * @param jnt_right
 */
void Robot_KIN::setCurrentJnt(JntArray jnt_left, JntArray jnt_right)
{
    jnt_current_left = jnt_left;
    jnt_current_right = jnt_right;
}

 /**
  * 将Frame中的数据拷贝到pose中
  * @param frame_in
  * @param pose_out
  */
void frame2pose(Frame frame_in, geometry_msgs::Pose &pose_out)
{
    pose_out.position.x = frame_in.p.x();
    pose_out.position.y = frame_in.p.y();
    pose_out.position.z = frame_in.p.z();

    frame_in.M.GetQuaternion(pose_out.orientation.x, pose_out.orientation.y, pose_out.orientation.z, pose_out.orientation.w);
}

 /**
  * 数据转换，将pose中的数据拷贝到Frame中去
  * @param pose_in      输入的pose
  * @param frame_out    输出的frame
  */
void pose2frame(geometry_msgs::Pose pose_in, Frame &frame_out)
{
    frame_out.p = Vector(pose_in.position.x, pose_in.position.y, pose_in.position.z);
    frame_out.M = Rotation::Quaternion(pose_in.orientation.x, pose_in.orientation.y
            , pose_in.orientation.z, pose_in.orientation.w);
}

 /**
  * 数据类型转换，将JntArray中的数据拷贝到std::vector<double>中去
  * @param jntarray_in      输入的JntArray
  * @param jntvector_out    输出的vector
  */
void jntArray2jntVector(JntArray jntarray_in, std::vector<double> &jntvector_out)
{
    for(unsigned int i = 0; i < jntvector_out.size(); i++)
    {
        jntvector_out[i] = jntarray_in.data[i];
    }
}

 /**
  * 将std::vector<double>中的数据拷贝到JntArray中去
  * @param jntvector_in     输入的vector
  * @param jntArray_out     输出的JntArray
  */
void jntVector2jntArray(std::vector<double> jntvector_in, JntArray &jntArray_out)
{
    for(unsigned int i = 0; i < jntvector_in.size(); i++)
    {
        jntArray_out.data[i] = jntvector_in[i];
    }
}

/**
 * 解析IK基础模块
 * @param frame_base    ： base连杆在world坐标系中的表示
 * @param frame_obj     ： 目标坐标系
 * @param q_sol         ： 得到的关节逆解
 * @param q6
 * @return              ： 返回解的数量
 */
inline int IK_analytical(Frame frame_base, Frame frame_obj, double *q_sol, double q6)
{
    double T[16] = {0};

    // ee_link到tool0的变换矩阵
    Frame frame_tool = Frame(Rotation(0, 0, 1, -1, 0, 0, 0, -1, 0));

    frame_obj = frame_base.Inverse() * frame_obj * frame_tool.Inverse();    // 求出极坐标系

    for(int i = 0; i < 3; i++)
    {
        // 获取旋转矩阵
        *(T + i) = frame_obj.M.data[i];
        *(T + 4 + i) = frame_obj.M.data[3 + i];
        *(T + 8 + i) = frame_obj.M.data[6 + i];

        *(T + 4*(i + 1) - 1) = frame_obj.p.data[i];     // 获取位置量
    }
    *(T + 15) = 1;

    return inverse(T, q_sol, q6);
}

/**
 * 选IK解。由于UR机器人基础解有8组，加上各关节运行范围为[-2*PI, 2*PI]，因此在要根据情况选择。这里主要根据离上一位置关节距离最近选择
 * TODO：三关节跳转处理
 * @param jnt_init  ： 初始关节角，选解时选择与此角度距离最近的解
 * @param q_sol     ： 逆解
 * @param q_num     ： 解的数量
 * @param jnt_out   ： 输出的关节变量
 * @return          ： 返回值，如果求解成功，则返回1，不成功返回-1
 */
inline int IK_selection(JntArray jnt_init, double *q_sol, int q_num, JntArray &jnt_out)
{
    double q_diff[8] = {0};       // 逆解到初始解的距离，此值辅助进行选解
    double min_dis = 100;   // 保存最小的关节距离值，辅助进行选解
    int min_index = 0;          // 最小距离的索引号，辅助进行选解
    double diff_temp;
  //  double q_w[6] = {4, 1, 1, 1, 1, 1};

 /*   double q_init[6];
    for(int i = 0; i < 6; i++)
        q_init[i] = jnt_init.data[i];
*/
    if(q_num > 0)
    {
        for (int i = 0; i < q_num; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                diff_temp = q_sol[i * 6 + j] - jnt_init.data[j];
                if(j != 2)  // 第三关节由于可能会导致碰撞，因此不能进行翻转
                {
                    if(diff_temp > PI && q_sol[i * 6 + j] > 0)  // 是否满足
                    {
                        q_sol[i * 6 + j] -= 2 * PI;
                        diff_temp -= 2.0 * PI;
                    }
                    else if(diff_temp < -PI && q_sol[i * 6 + j] < 0)
                    {
                        q_sol[i * 6 + j] += 2 * PI;
                        diff_temp += 2.0 * PI;
                    }
                }
                q_diff[i] += absv(diff_temp); // 计算各组逆解到关节初始值的距离
            }

            // 找到距离最小的那组解
            if (q_diff[i] < min_dis)
            {
                min_dis = q_diff[i];
                min_index = i;
            }
        }

        // 将距离最小的解保存到jnt_out中去
        for(int i = 0; i < 6; i++)
            jnt_out.data[i] = q_sol[min_index * 6 + i];

     //   if(q_diff[min_index] > 2*PI)
     //       return -1;

        if(absv(q_sol[min_index * 6] - jnt_init.data[0])> PI / 2)
 //           int k = 0;

        return 1;   // 求解成功，返回1
    } else {
        return -1;  // 求解失败，返回0
    }
}

/**
 * 初始化JointState变量，用来发布JointState消息，驱动Rviz机器人运动
 * @param initjnt ： 初始化值
 * @return ： 生成的JointState变量
 */
JointState InitJointState(JntArray initjnt)
{
    // 初始化
    JointState joint = JointState();
    joint.header = std_msgs::Header();
    joint.header.stamp = ros::Time::now();

    // 设定各关节名称
    joint.name.push_back("left_shoulder_pan_joint");
    joint.name.push_back("left_shoulder_lift_joint");
    joint.name.push_back("left_elbow_joint");
    joint.name.push_back("left_wrist_1_joint");
    joint.name.push_back("left_wrist_2_joint");
    joint.name.push_back("left_wrist_3_joint");

    joint.name.push_back("left_finger_joint_1");
    joint.name.push_back("left_finger_joint_2");
    joint.name.push_back("left_finger_joint_3");
//    joint.name.push_back("left_gripper_left_finger_joint");
//    joint.name.push_back("left_gripper_right_finger_joint");

    joint.name.push_back("right_shoulder_pan_joint");
    joint.name.push_back("right_shoulder_lift_joint");
    joint.name.push_back("right_elbow_joint");
    joint.name.push_back("right_wrist_1_joint");
    joint.name.push_back("right_wrist_2_joint");
    joint.name.push_back("right_wrist_3_joint");

    joint.name.push_back("right_finger_joint_2_1");
    joint.name.push_back("right_finger_joint_2_2");
//    joint.name.push_back("right_gripper_left_finger_joint");
//    joint.name.push_back("right_gripper_right_finger_joint");

    // 初始化关节角度
    joint.position.resize(joint.name.size());
    joint.velocity.resize(joint.name.size());

    for(unsigned int i = 0; i < joint.name.size(); i++)
    {

        joint.position[i] = initjnt.data[i];
    }
    return joint;
}

/**
 * 设置JointState中的关节值
 * @param jnt_left  ：左臂关节值
 * @param jnt_right ：右臂关节值
 * @param joint     ：JointState变量
 */
void setJointState(JntArray jnt_left, JntArray jnt_right, JointState &joint)
{
    for(int i = 0; i < 6; i++)
    {
        joint.position[i] = jnt_left.data[i];
        joint.position[i + 8] = jnt_right.data[i];
    }
}

/**
 * 设置JointState中的关节位置和速度
 * @param jnt_pos_left  ： 左臂关节位置
 * @param jnt_pos_right ： 右臂关节位置
 * @param jnt_vel_left  ： 左臂关节速度
 * @param jnt_vel_right ： 右臂关节速度
 * @param joint         ： JointState变量，输出值
 */
void setJointState(JntArray jnt_pos_left, JntArray jnt_pos_right, JntArray jnt_vel_left, JntArray jnt_vel_right, JointState &joint)
{
    for(int i = 0; i < 6; i++)
    {
        joint.position[i] = jnt_pos_left.data[i];
        joint.position[i + 8] = jnt_pos_right.data[i];

        joint.velocity[i] = jnt_vel_left.data[i];
        joint.velocity[i + 8] = jnt_vel_right.data[i];
    }
}

/**
 * 将Frame里的数据按照要求转换为一个float类型的数组，以便后续计算GMM的likelihood值
 * @param frame_arm     ： 输入的手臂末端的frame
 * @param data          ： 输出的5维的float数组
 */
void copyDataFromArmFrame(Frame frame_arm, float *data)
{
    *(data + 0) = (float)frame_arm.p.x();
    *(data + 1) = (float)frame_arm.p.y();
    *(data + 2) = (float)frame_arm.p.z();
    *(data + 3) = (float)frame_arm.M.UnitZ().x();
    *(data + 4) = (float)frame_arm.M.UnitZ().y();
}


/**
 * 由一个采样数据点获取得到中心坐标系
 * @param data          ： 采样数据点
 * @param frame_center  ： 输出得到的frame
 */
void getCenterFrameFromData(float *data, Frame & frame_center)
{
    Vector center_vector = Vector(*(data + 0), *(data + 1), *(data + 2));    // 中心点坐标

    Vector ref_vector = Vector(sqrt(3.0) / 3.0, sqrt(3.0) / 3.0, sqrt(3.0) / 3.0);
    Vector ref_vector1 = Vector(1.0, 0.0, 0.0);    // 参考的一个轴，用z轴和此轴求另外的几个轴坐标
    Vector Temp;

    float theta2 = *(data + 4);
    float phi2 = *(data + 5);

    Vector x_vector = Vector();         //
    Vector y_vector = Vector();
    Vector z_vector = Vector(sin(theta2) * cos(phi2), sin(theta2) * sin(phi2), cos(theta2));

    // 求x方向向量
    Temp = z_vector * ref_vector;
    if(Temp.Norm() < 10e-3)
    {
        Temp = z_vector * ref_vector1;
    }
    x_vector = Temp / Temp.Norm();
    y_vector = z_vector * x_vector;

    frame_center.M = Rotation(x_vector, y_vector, z_vector);
    frame_center.p = center_vector;

}

/**
 * 根据采样得到的六维数据，计算出各臂的位置
 * @param data          ： 六维位置数据（x, y, z, r, theta, phi）
 * @param frame_right   : 左臂末端点
 * @param frame_left    ：右臂末端点
 */
void getArmFrameFromData(float *data, Frame & frame_left, Frame & frame_right)
{
    Frame frame_center;

    float r2 = *(data + 3);
    getCenterFrameFromData(data, frame_center);
    getArmFrameFromCenter(frame_center, r2, frame_left, frame_right);
}

/**
 * 根据采样得到的六维向量，计算出各臂的位置
 * @param vector_in ：输入向量
 * @param frame_left ：输出左臂frame
 * @param frame_right ：输出右臂frame
 */
void getArmFrameFromVector(VectorXd vector_in, Frame &frame_left, Frame &frame_right)
{
    int dim = int(vector_in.rows());
    float *data = new float(dim);

    vectorxd2random(vector_in, data, dim);
  //  ROS_INFO("%f, %f, %f, %f, %f, %f", data[0], data[1], data[2], data[3], data[4], data[5]);
    getArmFrameFromData(data, frame_left, frame_right);
    delete data;
}

/**
 * 根据中心坐标位置，计算两手末端frame
 * @param frame_center  : 中心点Frame
 * @param r             ：两手距离
 * @param frame_left    ：左手末端Frame
 * @param frame_right   : 右手末端frame
 */
void getArmFrameFromCenter(Frame frame_center, float r, Frame &frame_left, Frame &frame_right)
{
    Vector v_x, v_y, v_z;

    v_x = frame_center.M.UnitX();
    v_y = frame_center.M.UnitY();
    v_z = frame_center.M.UnitZ();

    frame_left.p = frame_center.p + r * v_z;
    frame_left.M = Rotation(v_x, -v_y, -v_z);

    frame_right.p = frame_center.p - r * v_z;
    frame_right.M = Rotation(v_x, v_y, v_z);
}

/**
 * 根据中心位置、速度计算两臂末端的速度
 * @param frame_center ： 中心点
 * @param r            : 手臂到中心点距离
 * @param vel_center   : 中心点速度
 * @param vel_left     : 左臂末端速度
 * @param vel_right    : 右臂末端速度
 */
void getVelFromCenter(Frame frame_center, float r, Twist vel_center, Twist &vel_left, Twist &vel_right)
{
    Vector vector_left = -r * frame_center.M.UnitZ();
    Vector vector_right = r * frame_center.M.UnitZ();

    vel_left.rot = vel_center.rot;      // 两臂旋转速度与中心坐标系相同
    vel_right.rot = vel_center.rot;

    // 线速度 vm = vn + Pn,m * w;
    vel_left.vel = vel_center.vel + vector_left * vel_center.rot;
    vel_right.vel = vel_center.vel + vector_right * vel_center.rot;
}

void random2vectorxd(float *random_pos_in, Eigen::VectorXd & vector_out, int dim)
{
    for(int i = 0; i < dim; i++)
        vector_out[i] = *(random_pos_in + i);
}

void vectorxd2random(VectorXd vector_in, float *random_pos_out, int dim)
{
    for(int i = 0; i < dim; i++)
        *(random_pos_out + i) = vector_in[i];
}

void getCenterFrameFromVector(VectorXd vector, Frame & frame_center)
{
    float random[6];
    vectorxd2random(vector, random, 6);
    getCenterFrameFromData(random, frame_center);
}


Robot_Kinematics_Annalytical::Robot_Kinematics_Annalytical(KDL::Frame baseFrame_left, KDL::Frame baseFrame_right,
                                                           KDL::Frame toolFrame_left, KDL::Frame toolFrame_right)
{
    baseFrame_left_ = baseFrame_left;
    baseFrame_right_ = baseFrame_right;
    toolFrame_left_ = toolFrame_left;
    toolFrame_right_ = toolFrame_right;
}

/**
 * 无参数的构造函数，针对现在机器人的结构，可以简化一点
 */
Robot_Kinematics_Annalytical::Robot_Kinematics_Annalytical()
{
    baseFrame_left_ = KDL::Frame(KDL::Rotation::RPY(2.356194, 0, 3.1415926), KDL::Vector(0.0005, 0.22663, 1.53885));
    baseFrame_right_ = KDL::Frame(KDL::Rotation::RPY(2.356194, 0, 0), KDL::Vector(0.0005, -0.22663, 1.53885));
    toolFrame_left_ = KDL::Frame();
    toolFrame_right_ = KDL::Frame();
}

/**
 * 计算左臂的正向运动学
 * @param jnt_in
 * @param frame_out
 */
void Robot_Kinematics_Annalytical::FK_left(KDL::JntArray jnt_in, KDL::Frame &frame_out)
{
    double T[16] = {0};
    double q[6];

    Frame frame_ee_tool = Frame(Rotation(0, 0, 1, -1, 0, 0, 0, -1, 0));   // ee_link到tool0的变换矩阵

    KDL::Frame frame_temp = Frame();
    copyDataFromJntArrayToArray(jnt_in, q);
    forward(q, T);
    copyDataFromArrayToFrame(T, frame_temp);
    frame_out = baseFrame_left_ * frame_temp * frame_ee_tool * toolFrame_left_;
}

void Robot_Kinematics_Annalytical::FK_right(KDL::JntArray jnt_in, KDL::Frame &frame_out)
{
    double T[16] = {0};
    double q[6];
    Frame frame_ee_tool = Frame(Rotation(0, 0, 1, -1, 0, 0, 0, -1, 0));   // ee_link到tool0的变换矩阵

    KDL::Frame frame_temp = Frame();
    copyDataFromJntArrayToArray(jnt_in, q);
    forward(q, T);
    copyDataFromArrayToFrame(T, frame_temp);
    frame_out = baseFrame_right_ * frame_temp * frame_ee_tool * toolFrame_right_;
}

void Robot_Kinematics_Annalytical::IK_analytical_left(JntArray jnt_init, Frame frame, JntArray &jnt_out)
{
    double T[16] = {0};
    double q[8 * 6] = {0};
    int qnum;
    Frame frame_temp;

    Frame frame_ee_tool = Frame(Rotation(0, 0, 1, -1, 0, 0, 0, -1, 0));   // ee_link到tool0的变换矩阵
    frame_temp = baseFrame_left_.Inverse() * frame * frame_ee_tool.Inverse() * toolFrame_left_.Inverse();

    copyDataFromFrameToArray(frame_temp, T);

    qnum = inverse(T, q);

    IKSelection(q, qnum, jnt_init, jnt_out);
}

void Robot_Kinematics_Annalytical::IK_analytical_right(JntArray jnt_init, Frame frame, JntArray &jnt_out)
{
    double T[16] = {0};
    double q[8 * 6] = {0};
    int qnum;
    Frame frame_temp;

    Frame frame_ee_tool = Frame(Rotation(0, 0, 1, -1, 0, 0, 0, -1, 0));   // ee_link到tool0的变换矩阵
    frame_temp = baseFrame_right_.Inverse() * frame * frame_ee_tool.Inverse() * toolFrame_right_.Inverse();

    copyDataFromFrameToArray(frame_temp, T);

    qnum = inverse(T, q);

    IKSelection(q, qnum, jnt_init, jnt_out);
}

int Robot_Kinematics_Annalytical::IKSelection(double *q_sol, int q_num, KDL::JntArray jnt_init, KDL::JntArray &jnt_out)
{
    double q_diff[8] = {0};       // 逆解到初始解的距离，此值辅助进行选解
    double min_dis = 100;   // 保存最小的关节距离值，辅助进行选解
    int min_index = 0;          // 最小距离的索引号，辅助进行选解
    double diff_temp;

    if(q_num > 0)
    {
        for (int i = 0; i < q_num; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                diff_temp = q_sol[i * 6 + j] - jnt_init.data[j];
                if(j != 2)  // 第三关节由于可能会导致碰撞，因此不能进行翻转
                {
                    if(diff_temp > PI && q_sol[i * 6 + j] > 0)  // 是否满足
                    {
                        q_sol[i * 6 + j] -= 2 * PI;
                        diff_temp -= 2.0 * PI;
                    }
                    else if(diff_temp < -PI && q_sol[i * 6 + j] < 0)
                    {
                        q_sol[i * 6 + j] += 2 * PI;
                        diff_temp += 2.0 * PI;
                    }
                }
                q_diff[i] += absv(diff_temp); // 计算各组逆解到关节初始值的距离
            }

            // 找到距离最小的那组解
            if (q_diff[i] < min_dis)
            {
                min_dis = q_diff[i];
                min_index = i;
            }
        }

        // 将距离最小的解保存到jnt_out中去
        for(int i = 0; i < 6; i++)
            jnt_out.data[i] = q_sol[min_index * 6 + i];

        //   if(q_diff[min_index] > 2*PI)
        //       return -1;

       // if(absv(q_sol[min_index * 6] - jnt_init.data[0])> PI / 2)
            //           int k = 0;

        return 1;   // 求解成功，返回1
    }
    else
    {
        return -1;  // 求解失败，返回0
    }
}