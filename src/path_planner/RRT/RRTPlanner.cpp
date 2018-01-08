//
// Created by zhxin on 17-8-16.
//

#include "RRTPlanner.h"


/**
 * 采用RRT算法进行路径规划
 * @return
 */
PathSearch_result RRTPlanner::plan()
{
    Path_Node pathnode_temp(6);    // 采样用的path node

    // 创建一个保存连接树的向量
    VectorXd vector_out;
    vector_out.resize(6);

    Connection_result connection_result = connectCSpaceRRT(tree_start, _goal, vector_out);

    // 如果中间无碰撞，则可以直接退出了
    if(connection_result == connection_success)
    {
        path.push_back(_str);
        path.push_back(_goal);

        return pathsearch_finish;
    }

    int extendnum = 0;
    ExtendTree_result extend_result = extend_start_fail;

    // 一直扩展树，直到到达预定次数，或者规划成功
    while(extendnum < MAXEXTENDTIMES && extend_result != extend_finish)
    {
        // 扩展树
        extend_result = extendCSpaceTreesRRT(tree_start, tree_goal);
        ROS_INFO("%d", extendnum);
        extendnum++;
    }

    // 如果到达了规定的次数还没有规划成功，则返回规划失败
    if(extendnum == MAXEXTENDTIMES)
        return pathsearch_fail;

    std::cout << "采样点总数为：" << tree_start.size() + tree_goal.size() << std::endl;

    rearrangePath();

    std::cout << "压缩路径点前数量为： " << path.size() << std::endl;

    simplifyPath();

    std::cout << "压缩后数量为：" << path.size() << std::endl;

    return pathsearch_finish;
}


/**
 * 初始化类
 * @param str   ： 规划起始位置
 * @param goal  ： 规划目标位置
 */
void RRTPlanner::Init(VectorXd str, VectorXd goal)
{
    this->_str = str;
    this->_goal = goal;

    step = 0.01; // 碰撞检测的步长

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

    srand((unsigned)time(NULL));
}

/**
 * 给树添加节点操作
 * @param trees     ： 指定树结构
 * @param parent    ： 父节点编号
 * @param v         ： 节点处的向量
 */
void RRTPlanner::addTreeNode(std::vector<Path_Node> &trees, int parent, VectorXd v)
{
    Path_Node path_node_new;

    path_node_new.vector = v;
    path_node_new.parent = parent;

    VectorXd v_temp = v - trees[parent].vector;
    path_node_new.dis = v_temp.norm() + trees[parent].dis;

    int new_node_index = (int)trees.size();//新加的节点在树中的位置
    trees[parent].child.push_back(new_node_index);//更新父节点的子节点列表

//    ROS_INFO("The %d th node", (int)trees.size());
    trees.push_back(path_node_new);
}

/**
 * C空间进行连接到树操作
 * @param tree ：要连接的树
 * @param goal ：要连接的目标位置
 * @param vector_out ： 输出的向量位置
 * @return  ： 连接结果，如果连接成功，返回connection_success，如果一步都不能连接，则返回connection_fail，
 * 否则返回connection_mid，并将中间点的位置保存到vector_out中
 */
Connection_result RRTPlanner::connectCSpaceRRT(std::vector<Path_Node> &tree,
                                VectorXd goal, VectorXd &vector_out)
{
    int nearestIndex = findNearestNode(goal, tree);

    VectorXd str = tree[nearestIndex].vector;

    double len = 0;
    vector_out.setZero();

    // 获取单位向量
    VectorXd vector_temp = goal - str;
    VectorXd vector_dis_temp;
    VectorXd vectorUnitDirection = vector_temp / vector_temp.norm();

    bool iscollision = false;

    int addnum = 0;

    int parent_index;

    // 如果距离小于步长，则直接添加到树上
    if(vector_temp.norm() < step)
    {
        addTreeNode(tree, nearestIndex, goal);
        return connection_success;
    }

    // 以最小步长进行尝试连接，直到连接到了目标点，或者发生碰撞
    while((!iscollision) && len < vector_temp.norm())
    {
        len = len + step;
        vector_out = str + len * vectorUnitDirection;
        iscollision = collisionCheckRight(vector_out);

        // 如果不碰撞，则添加到树上去
        if(!iscollision)
        {
            parent_index = (addnum == 0) ? nearestIndex : (int)tree.size() - 1;
//            tree[parent_index]
            addTreeNode(tree, parent_index, vector_out);
            addnum++;
        }
    }

    // 如果两节点不发生碰撞，则返回0
    if (vector_temp.norm() - len < step)
    {
//        ROS_INFO("connection is free");
        return connection_success;
    }
    else if(len <= step)
    {
//        ROS_INFO("connection is invalid");
        return connection_fail;
    }
    else
    {
//        ROS_INFO("A mid point is found");
        // 最后一次保存的是碰撞的点，这里减去增量，得到最后无碰撞的点
        vector_out = vector_out - step * vectorUnitDirection;
        return connection_mid;
    }
}

Connection_result RRTPlanner::connectCSpaceRRT(VectorXd str, VectorXd goal)
{
    double len = 0;

    VectorXd vector_out;

    // 获取单位向量
    VectorXd vector_temp = goal - str;
    VectorXd vector_dis_temp;
    VectorXd vectorUnitDirection = vector_temp / vector_temp.norm();

    bool iscollision = false;

    // 如果距离小于步长，则直接添加到树上
    if(vector_temp.norm() < step)
    {
        return connection_success;
    }

    // 以最小步长进行尝试连接，直到连接到了目标点，或者发生碰撞
    while((!iscollision) && len < vector_temp.norm())
    {
        len = len + step;
        vector_out = str + len * vectorUnitDirection;
        iscollision = collisionCheckRight(vector_out);
    }

    // 如果两节点不发生碰撞，则返回0
    if (vector_temp.norm() - len < step)
        return connection_success;
    else
        return connection_fail;
}

/**
 * C空间中扩展树结构
 * @param str_tree ： 起始点处的树
 * @param goal_tree ： 目标点处的树
 * @return ：返回扩展结果
 */
ExtendTree_result RRTPlanner::extendCSpaceTreesRRT(std::vector<Path_Node> &str_tree, std::vector<Path_Node> &goal_tree)
{
    VectorXd vector_random, vector_temp_str, vector_temp_goal, dis_temp;
    Connection_result connection_result;

    Path_Node path_node_new(6);

    bool collision_flag = true;

    vector_random.resize(6);
    vector_temp_str.resize(6);
    vector_temp_goal.resize(6);

    // 随机生成一个无碰撞的点
    while(collision_flag) {
        for (int i = 0; i < 6; i++)
            vector_random[i] = getrandom(randomRangeLow[i], randomRangeUp[i]);

        collision_flag = collisionCheckRight(vector_random);
    }

//    std::cout << "random sample vector is: " <<vector_random.transpose() << std::endl;
    // 尝试进行连接
    connection_result = connectCSpaceRRT(tree_start, vector_random, vector_temp_str);

//    ROS_INFO("Start Tree Nodes: %d", tree_start.size());

    if(connection_result == connection_fail)
    {
//        ROS_INFO("extend start tree failed");
        return extend_start_fail;
    } // 如果不能连接，则返回失败
    else if(connection_result == connection_success)
        vector_temp_str = vector_random;    // 如果连接成功，则将最终连接的点设为这个随机点

//    ROS_INFO("extend start tree success");

    connection_result = connectCSpaceRRT(tree_goal, vector_temp_str, vector_temp_goal);

//    ROS_INFO("Goal Tree Nodes: %d", tree_goal.size());
    ExtendTree_result extend_result;

    if(connection_result == connection_fail)
    {
//        ROS_INFO("extend end fail");
        return extend_end_fail;
    }     // 如果连接不成功，则返回失败
    else if (connection_result == connection_success)
    {
//        ROS_INFO("Path found");
        extend_result = extend_finish;  // 如果连接上了，则返回已经完成
        vector_temp_goal = vector_temp_str;
    }
    else
    {
//        ROS_INFO("extend success");
        extend_result = extend_success;
    } // 如果连接上了一部分，则返回部分成功

    return extend_result;
}

ExtendTree_result RRTPlanner::extendCSpaceTreesRRTSingle(std::vector<Path_Node> &str_tree, VectorXd goal)
{
    VectorXd vector_random, vector_temp_str, dis_temp;
    Connection_result connection_result;

    Path_Node path_node_new(6);

    bool collision_flag = true;

    vector_random.resize(6);
    vector_temp_str.resize(6);

    VectorXd temp = goal - str_tree[0].vector;

    double dis = temp.norm();
    double sampleDis;

    // 随机生成一个无碰撞的点
    while(collision_flag)
    {
        for (int i = 0; i < 6; i++)
            vector_random[i] = getrandom(randomRangeLow[i], randomRangeUp[i]);

        sampleDis = (vector_random - str_tree[0].vector).norm() + (goal - vector_random).norm();

        // 这里只在小距离范围内采样
        if(sampleDis > 1.6 * dis)
            continue;

        collision_flag = collisionCheckRight(vector_random);
    }

    // 尝试进行连接
    connection_result = connectCSpaceRRT(str_tree, vector_random, vector_temp_str);

    ExtendTree_result extendTreeResult;
    if(connection_result == connection_fail)
    {
//        ROS_INFO("extend start tree failed");
        extendTreeResult = extend_start_fail;
    } // 如果不能连接，则返回失败
    else if(connection_result == connection_success)
    {
        extendTreeResult = extend_finish;
    }    // 如果连接成功，则将最终连接的点设为这个随机点
    else
        extendTreeResult = extend_success;

    return extendTreeResult;
}

/**
 * 找到树tree上距离p最近的节点的编号
 * @param p
 * @param tree
 * @return
 */
int RRTPlanner::findNearestNode(VectorXd p, std::vector<Path_Node> tree)
{
    double dis = 10000;
    int index = 0;
    VectorXd temp;

    for(unsigned int i = 0; i < tree.size(); i++)
    {
        temp = tree[i].vector - p;
        if(temp.norm() < dis)
        {
            index = i;
            dis = temp.norm();
        }
    }

    return index;
}

double RRTPlanner::rearrangePath()
{
    double totalCost = 0;   // 计算路径总代价

    // 先清除path
    path.clear();

    // 将规划的树重新排序到输出的变量rrt_path中去，并返回规划成功
    std::vector<VectorXd>::iterator it;

    int i = (int)tree_start.size() - 1;

    totalCost = tree_start[i].dis;

    // 将tree_start逆向排列
    while(i != 0)
    {
        it = path.begin();
        path.insert(it, tree_start[i].vector);
        i = tree_start[i].parent;
    }
    it = path.begin();
    path.insert(it, tree_start[0].vector);

    // 逐个添加tree_goal上的节点
    i = (int)tree_goal.size() - 1;

    totalCost += tree_goal[i].dis;
    while(i !=0 )
    {
        path.push_back(tree_goal[i].vector);
        i = tree_goal[i].parent;
    }
    path.push_back(tree_goal[0].vector);
    return totalCost;
}

RRTPlanner::RRTPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit):Planner(a, robMoveit)
{

}

RRTPlanner::~RRTPlanner(){}

void RRTPlanner::simplifyPath()
{
    Connection_result connectionResult;

    std::vector<VectorXd> pathTemp;

    if(path.size() > 2)
    {
        for(int i = 0; i < path.size() - 2; i++)
        {
            for(int j = (int)path.size() - 1; j > i + 2; j--)
            {
                connectionResult = connectCSpaceRRT(path[i], path[j]);
                if(connectionResult == connection_success)
                {
                    for(int kk = 0; kk <= i; kk++)
                        pathTemp.push_back(path[kk]);

                    for(int kk = j; kk < (int)path.size(); kk++)
                        pathTemp.push_back(path[kk]);

                    path = pathTemp;

                    simplifyPath();
                    return;
                }
            }
        }
    }
}