//
// Created by zhxin on 17-8-23.
//

#include "GMMGuidedPlanner.hpp"

#include "../GMM/GMM.h"

/**
 * 比较函数，排序的时候需要用到
 * @param a
 * @param b
 * @return
 */
inline bool compareDisClass(DisClass a, DisClass b)
{
    return a.dis < b.dis;
}

/**
 * 将数组里的数拷贝到VectorXd中去，在读取高斯均值的时候需要用到
 * @param arr
 * @param vector
 */


/**
 * 构造函数，初始化均值和方差
 */
myGmm::myGmm()
{
    CSpaceGmm = new Gmm(nstates, dim);

    for(int i = 0; i < nstates; i++)
    {
        CSpaceGmm->SetMean(i, mu[i]);                    // 设置均值
        CSpaceGmm->SetCovariance(i, sigma[i], false);   // 设置方差
        CSpaceGmm->SetPrior(i, prio[i]);              // 设置前验概率
    }

    thr = threshold;
}

/**
 * 根据
 * @param p
 * @return
 */
bool myGmm::checkValidationRight(VectorXd p)
{
    float *arrp;
    arrp = new float [CSpaceGmm->dim];
    float pdfValue;

    vector2arr(p, arrp);

  //  std::cout << p.transpose() << std::endl;

    pdfValue = CSpaceGmm->Pdf(arrp);

    pdfValue = log(double(pdfValue));

    delete[](arrp);

    return pdfValue > thr;
}

VectorXd myGmm::MahalanobisDis(VectorXd p)
{
    VectorXd dis;   // 定义距离向量
    dis.resize(CSpaceGmm->nstates);

    VectorXd mu;    // 均值
    mu.resize(CSpaceGmm->dim);

    VectorXd temp;

    MatrixXd cov;   // 方差矩阵
    cov.resize(CSpaceGmm->dim, CSpaceGmm->dim);

    float covArr[36];
    float muArr[6];

    for(int k = 0; k < CSpaceGmm->nstates; k++)
    {
        CSpaceGmm->GetCovariance(k, covArr, false); // 获取方差矩阵
        CSpaceGmm->GetMean(k, muArr);   // 获取均值矩阵
        for(int i = 0; i < CSpaceGmm->dim; i++)
        {
            for(int j = 0; j < CSpaceGmm->dim; j++)
            {
               cov(i, j) = covArr[i * CSpaceGmm->dim + j];  // 将方差矩阵保存到矩阵变量中
            }
            mu[i] = muArr[i];   // 将均值保存为向量
        }
        temp = p - mu;
        dis[k] = sqrt(temp.transpose() * cov.inverse() * temp); // 计算马氏距离
    }
    return dis;
}

/**
 * 利用马氏距离，找到附近最近的gaussian
 * @param p
 * @return
 */
int myGmm::findNearestGaussian(VectorXd p)
{
    VectorXd dis = angleDis(p);

//    VectorXd dis = MahalanobisDis(p);

 //   std::cout << dis.transpose() << std::endl;
    float dismin = 1000;
    int disminindex = 0;

    for(int i = 0; i < CSpaceGmm->nstates; i++)
    {
        if(dis[i] < dismin)
        {
            dismin = dis[i];
            disminindex = i;
        }
    }
    return disminindex;
}

VectorXd myGmm::angleDis(VectorXd p)
{
    VectorXd dis;
    dis.resize(CSpaceGmm->nstates);

    double disTemp;
    float *mui;
    mui = new float [CSpaceGmm->dim];
    for(int i = 0; i < CSpaceGmm->nstates; i++)
    {
        CSpaceGmm->GetMean(i, mui);   // 获取均值矩阵
        disTemp = 0;
        for(int j = 0; j < CSpaceGmm->dim; j++)
            disTemp += (mui[j] - p[j]) * (mui[j] - p[j]);
        dis[i] = sqrt(disTemp);
    }
    delete[](mui);

    return dis;
}

/**
 * 初始化
 */
void GMMGuidedPlanner::Init(VectorXd str, VectorXd goal)
{
    this->_str = str;
    this->_goal = goal;

    step = 0.05; // 碰撞检测的步长

    // 关节采样范围
    randomRangeLow.resize(6);
    randomRangeUp.resize(6);

    randomRangeLow << -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI;
    randomRangeUp << M_PI, M_PI, M_PI, M_PI, M_PI, M_PI;

    // 初始化树结构
    tree_start.clear();
    tree_goal.clear();

    Path_Node pathNodeTemp;
    pathNodeTemp.setNode(str, 0, 0);

    tree_start.push_back(pathNodeTemp);

    pathNodeTemp.vector = goal;
    tree_goal.push_back(pathNodeTemp);


    E.resize(GMM.CSpaceGmm->nstates, GMM.CSpaceGmm->nstates);
    disMat.resize(GMM.CSpaceGmm->nstates, GMM.CSpaceGmm->nstates);

    vexNum = GMM.CSpaceGmm->nstates;
}

/**
 * 找到最近的节点
 * @param p
 * @return
 */
int GMMGuidedPlanner::findNearestGaussian(VectorXd p)
{
    return GMM.findNearestGaussian(p);
}


VectorXd GMMGuidedPlanner::getCenter(int i)
{
    float *mu;
    mu = new float[GMM.CSpaceGmm->dim];
    VectorXd centerVector;
    centerVector.resize(GMM.CSpaceGmm->dim);

    GMM.CSpaceGmm->GetMean(i, mu);
    arr2vector(mu, centerVector);

    delete[](mu);
    return centerVector;
}

void GMMGuidedPlanner::constructRoadMap()
{
    int neighbor;

    Connection_result connection_result;

    std::vector<DisClass> disVector;

    disMat.setConstant(INFINITY);   // 将距离值设为最大值
    edgeNum = 0;

    E.setZero();

    float *muCenter, *muJ;
    muCenter = new float [GMM.CSpaceGmm->dim];
    muJ = new float [GMM.CSpaceGmm->dim];

    float disTemp;

    for(int i = 0; i < vexNum; i++)
    {
        GMM.CSpaceGmm->GetMean(i, muCenter);
        disVector = findSortedNeighborsProb(i);
        for(int j = 1; j < 5; j++)  // 从1开始，不计算自身
        {
            neighbor = disVector[j].index;
            if(E(i, neighbor) != 1)
            {
                GMM.CSpaceGmm->GetMean(neighbor, muJ);
                E(i, neighbor) = 1;
                E(neighbor, i) = 1;
                disTemp = 0;
                for(int kk = 0; kk < GMM.CSpaceGmm->dim; kk++)
                    disTemp += (muCenter[kk] - muJ[kk]) * (muCenter[kk] - muJ[kk]);

                disMat(i, neighbor) = (float)sqrt(disTemp);
                disMat(neighbor, i) = disMat(i, neighbor);

                edgeNum++;
            }
//            neighbor = disVector[j].index;
//            connection_result = connectCenter(i, neighbor); //进行判断是否可以直接相连，如果无碰撞相连，则添加E边和dis距离
//            if(connection_result == connection_success)
//            {
//                E(i, neighbor) = 1;
//                E(neighbor, i) = 1;
//                disMat(i, neighbor) = disVector[j].dis;
//                disMat(neighbor, i) = disVector[j].dis;
//                edgeNum++;
//            }
        }
        disMat(i, i) = 0;
    }

    delete [] muCenter;
    delete [] muJ;
//    std::cout << E << std::endl;
//
//    std::cout << disMat << std::endl;
}

std::vector<DisClass> GMMGuidedPlanner::findSortedNeighbors(int centerIndex)
{
    float *muCenter;
    float *muTemp;
    float disTemp;

    std::vector<DisClass> disVector;

    muCenter = new float [GMM.CSpaceGmm->dim];
    muTemp = new float [GMM.CSpaceGmm->dim];

    disVector.clear(); // 将原来数据清除

    DisClass disclass;



    // 分别计算centerrIndex中心点到其他中心点的距离
    GMM.CSpaceGmm->GetMean(centerIndex, muCenter);

    for(int i = 0; i < GMM.CSpaceGmm->nstates; i++)
    {
        disTemp = 0;
        GMM.CSpaceGmm->GetMean(i, muTemp);

        for(int j = 0; j < GMM.CSpaceGmm->dim; j++)
        {
            disTemp += (muTemp[j] - muCenter[j]) * (muTemp[j] - muCenter[j]);
        }
        disclass.dis = disTemp; // 因为只是排序，因此没有必要计算实际距离，利用平方来表示大小进行排序
        disclass.index = i;

        disVector.push_back(disclass);
    }

    // 进行排序
    std::sort(disVector.begin(), disVector.end(), compareDisClass);

    delete[](muCenter);
    delete[](muTemp);

    return disVector;
}

std::vector<DisClass> GMMGuidedPlanner::findSortedNeighborsProb(int centerIndex)
{
    float *muCenter;

    std::vector<DisClass> disVector;

    muCenter = new float [GMM.CSpaceGmm->dim];

    disVector.clear(); // 将原来数据清除

    DisClass disclass;

    // 分别计算centerrIndex中心点到其他中心点的距离
    GMM.CSpaceGmm->GetMean(centerIndex, muCenter);
    for(int i = 0; i < GMM.CSpaceGmm->nstates; i++)
    {
        disclass.dis = GMM.CSpaceGmm->Pdf(muCenter, i);
        disclass.dis = disclass.dis * GMM.CSpaceGmm->GetPrior(i);
        disclass.index = i;
        disVector.push_back(disclass);
    }

    // 进行排序
    std::sort(disVector.begin(), disVector.end(), compareDisClass);

    delete [] muCenter;

    return disVector;
}

/**
 * 计算第i个和第j个中心的距离
 * @param i
 * @param j
 * @return
 */
double GMMGuidedPlanner::calCenterDis(int i, int j)
{
    double dis = 0;
    float *mui, *muj;
    mui = new float [GMM.CSpaceGmm->dim];
    muj = new float [GMM.CSpaceGmm->dim];

    GMM.CSpaceGmm->GetMean(i, mui);
    GMM.CSpaceGmm->GetMean(j, muj);

    for(int k = 0; k < GMM.CSpaceGmm->dim; k++)
    {
        dis += (mui[k] - muj[k]) * (mui[k] - muj[k]);
    }

    dis = sqrt(dis);

    delete [] mui;
    delete [] muj;

    return dis;
}

Connection_result GMMGuidedPlanner::connectCenter(int i, int j)
{
    VectorXd mui, muj;
    mui.resize(GMM.CSpaceGmm->dim);
    muj.resize(GMM.CSpaceGmm->dim);

    float *arr;
    arr = new float[GMM.CSpaceGmm->dim];
    GMM.CSpaceGmm->GetMean(i, arr);
    arr2vector(arr, mui);
    GMM.CSpaceGmm->GetMean(j, arr);
    arr2vector(arr, muj);

//    std::cout << "Mui :" << mui.transpose() << std::endl;
//    std::cout << "Muj :" << muj.transpose() << std::endl;

    float ds = 0.1;
    float s = 0;

    bool collision = true;
    VectorXd direction = muj - mui;
    VectorXd temp;

    float dis = (float)direction.norm();
    direction = direction / dis;

    while(collision && s < dis)
    {
        s += ds;
        temp = mui + s * direction;
        collision = GMM.checkValidationRight(temp);
    }

    delete[](arr);

    if(s >= dis)
        return connection_success;
    else
        return connection_fail;
}

/**
 * 建立从strIndex处出发的最佳路径
 * @param strIndex
 * @return
 */
VectorXi GMMGuidedPlanner::findPathDigkstra(int strIndex)
{
    VectorXi flag;      // 标示节点是否被搜索过
    flag.resize(vexNum);
    flag.setZero();

    VectorXi preV;      // 前驱节点，即找到的最优路径中，上一个节点的index值
    preV.resize(vexNum);
    preV.setZero();

    VectorXd dist;      // 距离值
    dist.resize(vexNum);
    dist.setZero();

    float min;
    int nearestIndex = 0;
    float temp;

    // 从距离矩阵中，取出
    for (int i = 0; i < vexNum; i++)
    {
        dist[i] = disMat(strIndex, i);

        // 设置前向节点
        if(dist[i] == INFINITY)
            preV[i] = -1;
        else
            preV[i] = strIndex;
    }

    flag[strIndex] = 1;
    dist[strIndex] = 0;

  //  std::cout << E << std::endl;
    //std::cout << "dist vector is:" << dist.transpose() << std::endl;
    // 找到最近的点
    for(int i = 0; i < vexNum; i++)
    {
   //     std::cout << "dist:" << dist.transpose() << std::endl;
        min = INFINITY;
        for(int j = 0; j < vexNum; j++)
        {
            if(flag[j] == 0 && dist[j] < min)
            {
                min = dist[j];
                nearestIndex = j;
            }
        }

        flag[nearestIndex] = 1;

        for(int j = 0; j < vexNum; j++)
        {
            if(flag[j] == 0 && disMat(nearestIndex, j) < INFINITY)
            {
                if(dist[nearestIndex] + disMat(nearestIndex, j) < dist[j])
                {
                    dist[j] = dist[nearestIndex] + disMat(nearestIndex, j);
                    preV[j] = nearestIndex;
                }
            }
        }
    }
    return preV;
}

PathSearch_result GMMGuidedPlanner::plan()
{
    int strIndex = findNearestGaussian(_str);
    int goalIndex = findNearestGaussian(_goal);

    std::vector<Path_Node> strTree, goalTree;   // 从起始点规划到第一个高斯中心，和目标点

    std::vector<VectorXd> plannedPath;

    PathSearch_result pathSearchResult;
    //第一步，先规划str到第strIndex中心点的路径，此路径采用标准的单向RRT算法
    pathSearchResult = singleRRTSearch(_str, strIndex, strTree);
    if(pathSearchResult == pathsearch_fail)
        return pathsearch_fail;

    //将第一段路径添加到已规划好的路径中
    std::vector<VectorXd>::iterator it;
    int insertIndex = (int)strTree.size() - 1;
    while(insertIndex != 0)
    {
        it = plannedPath.begin();
        plannedPath.insert(it, strTree[insertIndex].vector);
        insertIndex = strTree[insertIndex].parent;
    }
    it = plannedPath.begin();
    plannedPath.insert(it, strTree[0].vector);
    path.clear();
    path = plannedPath;
    simplifyPath();
    plannedPath = path;

    // 第二步，规划goal到第goalIndex中心点的路径，也采用标准的单向RRT算法
    pathSearchResult = singleRRTSearch(_goal, goalIndex, goalTree);
    if(pathSearchResult == pathsearch_fail)
        return pathsearch_fail;

    // 第三步，规划strIndex到goalIndex中心的路径，此部分采用双向RRT算法，并且在利用PRM规划出来的最优路径附近采样。
    VectorXi preV;
    preV = findPathDigkstra(strIndex);    // 找到最优的中心路径
 //   std::cout << preV << std::endl;
    std::vector<int> sampleIndexVector;

    int addvexIndex = goalIndex;
    sampleIndexVector.push_back(addvexIndex);
    while(addvexIndex != strIndex)
    {
        addvexIndex = preV[addvexIndex];
        sampleIndexVector.push_back(addvexIndex);
        std::cout << "GMM states:" << addvexIndex <<std::endl;
    }


    // 第四部，检验规划的路径是否存在碰撞，如果存在碰撞，则在相应区域重新规划
    float *mu;
    mu = new float[GMM.CSpaceGmm->dim];
    VectorXd centerVectorXd;
    centerVectorXd.resize(GMM.CSpaceGmm->dim);
    bool *centerCollisionFlag;
    centerCollisionFlag = new bool[sampleIndexVector.size()];     // 高斯球中心点是否碰撞标志位
    bool *edgeCollisionFlag;
    edgeCollisionFlag = new bool[sampleIndexVector.size() - 1];   // 各边是否碰撞标志

    // 首先检测各高斯球中心点是否存在碰撞
    for(int i = 0; i < sampleIndexVector.size(); i++)
    {
        GMM.CSpaceGmm->GetMean(sampleIndexVector[i], mu);
        arr2vector(mu, centerVectorXd);
        centerCollisionFlag[i] = collisionCheckRight(centerVectorXd);   // 保存第路径中第i个节点是否存在碰撞
    }

    double ds = 0.04, s = 0, dis;
    VectorXd strVector, goalVector, temp, unitVector;
    std::vector<int> planningCenters;
    std::vector<VectorXd> innerPath;

    bool collisionFlag = false;
    // 如果第i各点处发生了碰撞，那么检测i-1及i+1是否发生碰撞，如果没有，那么重新规划 i - 1 到 i + 1段
    for(int i = 1; i < sampleIndexVector.size() - 1; i++)
    {
        innerPath.clear();
        if(centerCollisionFlag[i])  // 如果发生碰撞，则重新规划i - 1 到 i + 1段
        {
            planningCenters.clear();
            planningCenters.push_back(i);

            // 找到碰撞的区段，对碰撞的上一点到下一点中间段进行规划，并保存相应的高斯编号，采样时着重在此区域采样
            while(centerCollisionFlag[i])
            {
                i++;
                planningCenters.push_back(i);
            }

            // 重新规划这一段
            planningInGMM(planningCenters, innerPath);

            // 将路径保存下来
            for(int k = 0; k < innerPath.size(); k++)
                plannedPath.push_back(innerPath[k]);

        }
        else
        {
            strVector = getCenter(i - 1);
            goalVector = getCenter(i);
            unitVector = goalVector - strVector;
            dis = unitVector.norm();
            unitVector = unitVector / dis;
            while(s < dis)
            {
                s = s + ds;
                temp = strVector + s * unitVector;
                collisionFlag = collisionCheckRight(temp);
                if(collisionFlag)
                    break;
            }

            if(collisionFlag)
            {
                planningCenters.push_back(i - 1);
                planningCenters.push_back(i);
                planningInGMM(planningCenters, innerPath);

                for(int k = 1; k < innerPath.size(); k++)
                    plannedPath.push_back(innerPath[k]);
            } else
            {
                plannedPath.push_back(goalVector);
            }

        }
    }

    path.clear();
    path = plannedPath;
    // 逐个添加tree_goal上的节点
    insertIndex = (int)goalTree.size() - 1;
    while(insertIndex !=0 )
    {
        path.push_back(goalTree[insertIndex].vector);
        insertIndex = goalTree[insertIndex].parent;
    }
    path.push_back(goalTree[0].vector);

    delete[](centerCollisionFlag);
    delete[](edgeCollisionFlag);
    delete[](mu);

    return pathsearch_finish;
}


PathSearch_result GMMGuidedPlanner::singleRRTSearch(VectorXd p, int centerIndex, std::vector<Path_Node> &tree)
{
    float *muCenter;
    muCenter = new float [GMM.CSpaceGmm->dim];

    VectorXd centerVector;
    centerVector.resize(GMM.CSpaceGmm->dim);

    // 获取起始点和终点的位置值
    GMM.CSpaceGmm->GetMean(centerIndex, muCenter);
    arr2vector(muCenter, centerVector);
    delete[](muCenter);

    // 添加初始节点到树上
    Path_Node pathNode;
    pathNode.vector = p;
    pathNode.parent = -1;   // -1表示此节点为初始节点
    pathNode.dis = 0;
    tree.push_back(pathNode);

    VectorXd temp = centerVector - p;
    float disStart2Goal = temp.norm();

    Connection_result connectionResult = connectCSpaceRRT(p, centerVector);
    // 如果可以直接连接，则说明可以直接规划过去，那么就不需要添加sampling点了
    if(connectionResult == connection_success)
    {
        addTreeNode(tree, 0, centerVector);
        return pathsearch_finish;
    }

    // 如果不能直接连接，则采用有限范围内的sampling extend tree
    ExtendTree_result extendTreeResult = extend_start_fail;
    int extendTimes = 0;    // 搜索次数
    int maxExtendTimes = 10000; // 最大搜索次数
    while(extendTreeResult != extend_finish && extendTimes < maxExtendTimes)
    {
        extendTreeResult = extendCSpaceTreesRRTSingle(tree, centerVector);
        extendTimes++;
    }

    // 如果次数到了还没成功，则返回路径规划失败
    if(extendTimes >= maxExtendTimes)
        return pathsearch_fail;
    else    // 否则是规划成功
        return pathsearch_finish;
}

PathSearch_result GMMGuidedPlanner::planningInGMM(std::vector<int> centers, std::vector<VectorXd> &path)
{
    std::vector<Path_Node> strTree, goalTree;
    Path_Node pathNode;

    // 获取初始位置和末端位置
    VectorXd strVector = getCenter(centers[0]);
    VectorXd goalVector = getCenter(centers[centers.size() - 1]);

    // 将起始点和终点添加到树结构中
    pathNode.vector = strVector;
    pathNode.dis = 0;
    pathNode.parent = -1;
    strTree.push_back(pathNode);

    pathNode.vector = goalVector;
    pathNode.dis = 0;
    pathNode.parent = -1;
    goalTree.push_back(pathNode);


    // 创建一个保存连接树的向量
    VectorXd vector_out;
    vector_out.resize(6);

    Connection_result connection_result = connectCSpaceRRT(strTree, goalVector, vector_out);

    // 如果中间无碰撞，则可以直接退出了
    if(connection_result == connection_success)
    {
        path.push_back(strVector);
        path.push_back(goalVector);
        return pathsearch_finish;
    }

    int extendnum = 0;
    ExtendTree_result extend_result = extend_start_fail;

    // 一直扩展树，直到到达预定次数，或者规划成功
    while(extendnum < MAXEXTENDTIMES && extend_result != extend_finish)
    {
        // 扩展树
        extend_result = extendCSpaceTreesRRTInGMM(strTree, goalTree, centers);
        extendnum++;
    }

    // 如果到达了规定的次数还没有规划成功，则返回规划失败
    if(extendnum == MAXEXTENDTIMES)
        return pathsearch_fail;

    rearrangePath();

    path = this->path;

    return pathsearch_finish;
}


ExtendTree_result GMMGuidedPlanner::extendCSpaceTreesRRTInGMM(std::vector<Path_Node> &str_tree
        , std::vector<Path_Node> &goal_tree
        , std::vector<int> centers)
{
    VectorXd vector_random, vector_temp_str, vector_temp_goal, dis_temp;
    Connection_result connection_result;

    Path_Node path_node_new(6);

    bool collision_flag = true;

    vector_random.resize(6);
    vector_temp_str.resize(6);
    vector_temp_goal.resize(6);

    float *randomArr;
    randomArr = new float [GMM.CSpaceGmm->dim];
    long randomInt;
    int randomCenter;
    // 随机生成一个无碰撞的点
    while(collision_flag)
    {
        // 在centers中随机选择一个高斯
        randomInt = random();
        randomCenter = centers[randomInt % centers.size()];
        GMM.CSpaceGmm->Pdf(randomArr, randomCenter);    // 再利用这个高斯随机生成一个点
        arr2vector(randomArr, vector_random);       // 转换到vector中
        collision_flag = collisionCheckRight(vector_random);    // 判断此点是否存在碰撞， 如果碰撞，那么继续生成
    }

    delete[](randomArr);

    // 尝试进行连接
    connection_result = connectCSpaceRRT(str_tree, vector_random, vector_temp_str);

    if(connection_result == connection_fail)
    {
        ROS_INFO("extend start tree failed");
        return extend_start_fail;
    } // 如果不能连接，则返回失败
    else if(connection_result == connection_success)
        vector_temp_str = vector_random;    // 如果连接成功，则将最终连接的点设为这个随机点

    ROS_INFO("extend start tree success");

    connection_result = connectCSpaceRRT(goal_tree, vector_temp_str, vector_temp_goal);

    ExtendTree_result extend_result;

    if(connection_result == connection_fail)
    {
        ROS_INFO("extend end fail");
        return extend_end_fail;
    }     // 如果连接不成功，则返回失败
    else if (connection_result == connection_success)
    {
        ROS_INFO("Path found");
        extend_result = extend_finish;  // 如果连接上了，则返回已经完成
        vector_temp_goal = vector_temp_str;
    }
    else
    {
        ROS_INFO("extend success");
        extend_result = extend_success;
    } // 如果连接上了一部分，则返回部分成功

    return extend_result;
}


void update()
{

}