//
// Created by zhxin on 16-12-19.
//
#include "myGmm.h"
#include "../GMM/GMM_left.h"
#include "../GMM/GMM_right.h"
#include "../GMM/GMM.h"
/**
 * 由于训练的GMM的数据是经过归一化的，那么测试的数据也要进行相应计算才能求likelihood
 * @param data_in          : 输入的数据
 * @param preprocessing    ： 最大最小值，归一化需要用到的值
 * @param randcol          ： 输出数据
 * @param dim              : 数据维度
 * @param data_out         ： 输出的数据
 */
void normalizedata(float *data_in, float *preprocessing, int *randcol, int dim, float *data_out)
{

    float temp[DIMENSION] = {0};
    for(int i = 0; i < dim; i++)
    {
        temp[i] = (*(data_in + i) - *(preprocessing + i)) / (*(preprocessing + dim + i) - *(preprocessing + i));
    }

    for(int i = 0; i < dim; i++)
    {
        *(data_out + i) = temp[randcol[i] - 1];
    }
}


/**
 * 由于训练的GMM的数据是经过归一化的，那么测试的数据也要进行相应计算才能求likelihood
 * @param vector_in : 输入的数据
 * @param preprocessing ：最大最小值
 * @param randcol ： 随机列号
 * @param dim   ： 维度
 * @param vector_out ：输出的向量
 */
void normalizedata(VectorXd vector_in, float *preprocessing, int *randcol, int dim, VectorXd &vector_out)
{
    float *data_in = new float(dim);
    float *data_out = new float(dim);
    vectorxd2random(vector_in, data_in, dim);
    normalizedata(data_in, preprocessing, randcol, dim, data_out);
    random2vectorxd(data_out, vector_out, dim);
    delete data_in;
    delete data_out;
}


/**
 * 由一个采样数据点获取得到中心坐标系
 * @param data          ： 采样数据点
 * @param frame_center  ： 输出得到的frame
 */
/*void getCenterFrameFromData(float *data, Frame & frame_center)
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
*/
/**
 * 根据采样得到的六维数据，计算出各臂的位置
 * @param data          ： 六维位置数据（x, y, z, r, theta, phi）
 * @param frame_right   : 左臂末端点
 * @param frame_left    ：右臂末端点
 */
/*void getArmFrameFromData(float *data, Frame & frame_left, Frame & frame_right)
{
    Frame frame_center;

    float r2 = *(data + 3);
    getCenterFrameFromData(data, frame_center);
    getArmFrameFromCenter(frame_center, r2, frame_left, frame_right);
}
*/
/**
 * 根据中心坐标位置，计算两手末端frame
 * @param frame_center  : 中心点Frame
 * @param r             ：两手距离
 * @param frame_left    ：左手末端Frame
 * @param frame_right   : 右手末端frame
 */
/*void getArmFrameFromCenter(Frame frame_center, float r, Frame &frame_left, Frame &frame_right)
{
    Vector v_x, v_y, v_z;

    v_x = frame_center.M.UnitX();
    v_y = frame_center.M.UnitY();
    v_z = frame_center.M.UnitZ();

    frame_left.p = frame_center.p + r * v_z;
    frame_left.M = Rotation(v_x, -v_y, -v_z);

    frame_right.p = frame_center.p - r * v_z;
    frame_right.M = Rotation(v_x, v_y, v_z);
}*/

/**
 * 将Frame里的数据按照要求转换为一个float类型的数组，以便后续计算GMM的likelihood值
 * @param frame_arm     ： 输入的手臂末端的frame
 * @param data          ： 输出的5维的float数组
 */
 /*
inline void copyDataFromArmFrame(Frame frame_arm, float *data)
{
    *(data + 0) = (float)frame_arm.p.x();
    *(data + 1) = (float)frame_arm.p.y();
    *(data + 2) = (float)frame_arm.p.z();
    *(data + 3) = (float)frame_arm.M.UnitZ().x();
    *(data + 4) = (float)frame_arm.M.UnitZ().y();
}*/

/**
 * 根据中心位置、速度计算两臂末端的速度
 * @param frame_center ： 中心点
 * @param r            : 手臂到中心点距离
 * @param vel_center   : 中心点速度
 * @param vel_left     : 左臂末端速度
 * @param vel_right    : 右臂末端速度
 */
 /*
void getVelFromCenter(Frame frame_center, float r, Twist vel_center, Twist &vel_left, Twist &vel_right)
{
    Vector vector_left = -r * frame_center.M.UnitZ();
    Vector vector_right = r * frame_center.M.UnitZ();

    vel_left.rot = vel_center.rot;      // 两臂旋转速度与中心坐标系相同
    vel_right.rot = vel_center.rot;

    // 线速度 vm = vn + Pn,m * w;
    vel_left.vel = vel_center.vel + vector_left * vel_center.rot;
    vel_right.vel = vel_center.vel + vector_right * vel_center.rot;
}*/

Robot_Gmm::Robot_Gmm()
{
    // 初始化GMM模型
    Gmm_dual = new Gmm(35, 6);
    Gmm_left = new Gmm(30, 5);
    Gmm_right = new Gmm(30, 5);

    for(int i = 0; i < NBSTATES; i++) {
        Gmm_dual->SetMean(i, mu[i]);                    // 设置均值
        Gmm_dual->SetCovariance(i, sigma[i], false);   // 设置方差
        Gmm_dual->SetPrior(i, prio[i]);              // 设置前验概率
    }

    // 设置两臂各自的GMM参数
    for(int i = 0; i < 30; i++)
    {
        Gmm_left->SetMean(i, mu_left[i]);
        Gmm_left->SetPrior(i, prio_left[i]);
        Gmm_left->SetCovariance(i, sigma_left[i], false);

        Gmm_right->SetMean(i, mu_right[i]);
        Gmm_right->SetPrior(i, prio_right[i]);
        Gmm_right->SetCovariance(i, sigma_right[i], false);
    }

    // likelihood阈值
    thre =  4.6856;
    thre_left = 5.00;
    thre_right =  4.9067;
}

Robot_Gmm::~Robot_Gmm() {}

bool Robot_Gmm::isValidDualArm(VectorXd vector)
{
    VectorXd normvector;
    normvector.resize(6);
    normalizedata(vector, preprocessing[0], randcol, 6, normvector); // 数据归一化
    float likelihood = Gmm_dual->Pdf(normvector); // 计算此位置处的likelihood值，如果大于阈值，则认为此点时可行的
    return likelihood >= thre;
}

bool Robot_Gmm::isValidLeftArm(VectorXd vector)
{
    Frame frame_obj_left = Frame();
    Frame frame_obj_right = Frame();

    float normdata_left[5];
    getArmFrameFromVector(vector, frame_obj_left, frame_obj_right);

    // 计算左臂的可达性
    copyDataFromArmFrame(frame_obj_left, normdata_left);
    normalizedata(normdata_left, preprocessing_left[0], randcol_left, 5, normdata_left);
    float likelihood_left = Gmm_left->Pdf(normdata_left);
    return likelihood_left >= thre_left;
}

bool Robot_Gmm::isValidRightArm(VectorXd vector) {
    Frame frame_obj_left = Frame();
    Frame frame_obj_right = Frame();

    float normdata_right[5];
    getArmFrameFromVector(vector, frame_obj_left, frame_obj_right);

    // 计算左臂的可达性
    copyDataFromArmFrame(frame_obj_right, normdata_right);
    normalizedata(normdata_right, preprocessing_right[0], randcol_right, 5, normdata_right);
    float likelihood_right = Gmm_right->Pdf(normdata_right);
    return likelihood_right >= thre_right;
}

bool Robot_Gmm::isValidAll(VectorXd vector) {
    bool isvalid_dual, isvalid_left, isvalid_right;
    isvalid_dual = isValidDualArm(vector);
    isvalid_left = isValidLeftArm(vector);
    isvalid_right = isValidRightArm(vector);
    return isvalid_dual && isvalid_left && isvalid_right;
}