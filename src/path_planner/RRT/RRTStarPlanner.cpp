//
// Created by zhxin on 17-8-17.
//

#include "RRTStarPlanner.h"

RRTStarPlanner::RRTStarPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit):RRTPlanner(a, robMoveit)
{
    radius = 3.0;
}

RRTStarPlanner::~RRTStarPlanner(){}

PathSearch_result RRTStarPlanner::plan()
{
    Path_Node pathnode_temp(6);    // 采样用的path node
    std::vector<Feasible_Path> feasible_path;//保存所有可行路径

    // 创建一个保存连接树的向量
    VectorXd vector_out;
    vector_out.resize(6);

//    Connection_result connection_result = connectCSpaceRRT(tree_start, _goal, vector_out);
//
//    // 如果中间无碰撞，则可以直接退出了
//    if(connection_result == connection_success)
//    {
//        path.push_back(_str);
//        path.push_back(_goal);
//
//        return pathsearch_finish;
//    }

    int extendnum = 0;  //扩展次数
    ExtendTree_result extend_result = extend_start_fail;

    int feasible_path_num = 0;
    double total_cost;
    //寻找并保存可行路径
    Feasible_Path feasible;
    while(feasible_path_num < 8)
    {
        // 一直扩展树，直到到达预定次数，或者规划成功
        while(extendnum < MAXEXTENDTIMES && extend_result != extend_finish)
        {
            // 扩展树
            extend_result = extendCSpaceTreesRRTStar(tree_start, tree_goal);
            extendnum++;
        }
        // 保存可行路径及其cost
        if(extend_result == extend_finish)
        {
            total_cost = rearrangePath();
            feasible.saveFeasiblePath(path, total_cost);
            feasible_path.push_back(feasible);
            feasible_path_num++;
        }
        extendnum = 0;//清零
        extend_result = extend_start_fail;//重置
    }


//    // 如果到达了规定的次数还没有规划成功，则返回规划失败
//    if(extendnum == MAXEXTENDTIMES)
//        return pathsearch_fail;
//
//    rearrangePath();
//
//    simplifyPath();

    double totalCostMin = 10000.0;
    std::ofstream  outCostData;
    std::ofstream  outOptimalPathData;
    outCostData.open("/home/geds/catkin_rrt/src/dual_arm_robot/src/path_planning/data/cost.txt");
    outOptimalPathData.open("/home/geds/catkin_rrt/src/dual_arm_robot/src/path_planning/data/OptimalPath.txt");
    for(int i=0; i<feasible_path_num; i++)
    {
        if(feasible_path[i].total_cost < totalCostMin)
        {
            optimalPath.clear();//清零
            totalCostMin = feasible_path[i].total_cost;
            optimalPath = feasible_path[i].feasiblepath;
//            outCostData << optimalPath[i] << std::endl;
        }
        outCostData << feasible_path[i].total_cost << std::endl;
    }
    outCostData << "totalCostMin = " << totalCostMin << std::endl;
    int lenTemp = (int)optimalPath.size();
    for(int i = 0; i < lenTemp; i++)
    {
        outOptimalPathData << optimalPath[i].transpose() << std::endl;//保存最优路径
    }
    outCostData.close();
    outOptimalPathData.close();

    return pathsearch_finish;


    return pathsearch_finish;
}


ExtendTree_result RRTStarPlanner::extendCSpaceTreesRRTStar(std::vector<Path_Node> &str_tree, std::vector<Path_Node> &goal_tree)
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

    // 尝试进行连接
    connection_result = connectCSpaceRRT(tree_start, vector_random, vector_temp_str);

    if(connection_result == connection_fail)
    {
        return extend_start_fail;
    } // 如果不能连接，则返回失败
    else if(connection_result == connection_success)
    {
        vector_random = vector_temp_str;//改1
        //vector_temp_str = vector_random;    // 如果连接成功，则将最终连接的点设为这个随机点
    }


    //ChooseBestParent,选择总路径最短的父节点
    std::vector<int> neighbors_others;

    // 找到最佳的节点以及周边近的节点
    int bestParentIndex = chooseBestParent(vector_temp_str, tree_start, neighbors_others);

    addTreeNode(tree_start, bestParentIndex, vector_temp_str);  // 将最近点添加为父节点

    int newVectorIndex = (int)tree_start.size() - 1;    // 新添加的点的index值

    std::vector<int> rewireNeighbors;

    // 重新连接附近的节点
    rewirePathNodes(vector_temp_str, tree_start, neighbors_others, rewireNeighbors, newVectorIndex);

    vector_random = tree_start[rewireNeighbors[0]].vector;//将str_tree上rewire的节点向量（若有多个，只用第一个）赋给vector_random

    // 扩展goal端的树
    connection_result = connectCSpaceRRT(tree_goal, vector_random, vector_temp_goal);

    if(connection_result == connection_fail)
    {
        ROS_INFO("extend goal tree failed");
        return extend_end_fail;//本次扩展结束，重新扩展
    } // 如果不能连接，则返回失败，并重新生成随机点，重新连接
    else if(connection_result == connection_success)
    {
        ROS_INFO("Path found");
        return extend_finish;  // 如果连接上了，则返回已经完成
    }


    // 找到最佳的节点以及周边近的节点
    bestParentIndex = chooseBestParent(vector_temp_goal, tree_goal, neighbors_others);

    addTreeNode(tree_goal, bestParentIndex, vector_temp_goal);  // 将最近的点添加为父节点

    newVectorIndex = (int)tree_goal.size() - 1;

    rewireNeighbors.clear();

    //优化附近的其他节点
    rewirePathNodes(vector_temp_goal, tree_goal, neighbors_others, rewireNeighbors, newVectorIndex);

    return extend_success;
}

Connection_result RRTStarPlanner::connectCSpaceRRTStar(std::vector<Path_Node> &tree, VectorXd goal, VectorXd &vector_out, int nodeIndex)
{
//    int nearestIndex = findNearestNode(goal, tree);
    int nearestIndex = nodeIndex;
    VectorXd str = tree[nearestIndex].vector;

    double len = step;

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
//        len = len + step;
        vector_out = str + len * vectorUnitDirection;
        iscollision = collisionCheckRight(vector_out);

        // 如果不碰撞，则添加到树上去
        if(!iscollision)
        {
       //     parent_index = (addnum == 0) ? nearestIndex : (int)tree.size() - 1;
       //     addTreeNode(tree, parent_index, vector_out);
            len = len + step;
            addnum++;
        }
    }

    // 如果两节点不发生碰撞，则返回0
    if (vector_temp.norm() - len < step)
    {
        ROS_INFO("connection is free");
        return connection_success;
    }
    else if(len <= step)
    {
        ROS_INFO("connection is invalid");
        return connection_fail;
    }
    else
    {
        ROS_INFO("A mid point is found");
        // 最后一次保存的是碰撞的点，这里减去增量，得到最后无碰撞的点
        vector_out = vector_out - step * vectorUnitDirection;
        return connection_mid;
    }
}

void RRTStarPlanner::rewirePathNodes(VectorXd vector_new, std::vector<Path_Node> trees, std::vector<int> &neighbors_others,
                     std::vector<int> &rewire_neighbors, int vector_new_index)
{
    VectorXd temp_vector_differ;
    double dis_init2new2neighbor;//起始点到vector_new的距离
    int length = (int)neighbors_others.size();

    int k = 0;//能继续优化的邻居的个数
    int child_index_temp;//子节点索引号
    for(int i=0; i < length; i++)
    {
        temp_vector_differ = trees[neighbors_others[i]].vector - vector_new;
        dis_init2new2neighbor = temp_vector_differ.norm() + trees[vector_new_index].dis;
        if(dis_init2new2neighbor < trees[neighbors_others[i]].dis )
        {
            trees[neighbors_others[i]].parent = vector_new_index;//若距离优化了，将vector_new_index作为该neighbor的父节点索引号
            trees[neighbors_others[i]].dis = dis_init2new2neighbor;//同时更新距离
            //更新可继续优化的neighbor的所有子节点(如果有的话)的cost
            if(!trees[neighbors_others[i]].child.empty())
            {
                //              std::cout << trees[neighbors_others[i]].child.size() << std::endl;
                modifyChildCost(trees, neighbors_others[i], neighbors_others[i]);
            }
//            int len = (int)trees[neighbors_others[i]].child.size();
//            for(int j=0; j<len; j++)
//            {
//                child_index_temp = trees[neighbors_others[i]].child[j];
//                temp_vector_differ = trees[neighbors_others[i]].vector - trees[child_index_temp].vector;
//                trees[child_index_temp].dis = dis_init2new2neighbor + temp_vector_differ.norm();//更新子节点的cost
//            }
            rewire_neighbors.push_back(neighbors_others[i]);//保存能优化的neighbors在树中的位置
            k++;
        }
    }
    if(k==0)
        rewire_neighbors.push_back(vector_new_index);//如果rewire失败，则返回vector_new_index
}

int RRTStarPlanner::chooseBestParent(VectorXd vector_new, std::vector<Path_Node> trees, std::vector<int> &neighbors_others)
{
    neighbors_others.clear();

    VectorXd temp_differ;
    Connection_result connect_neighber_result;
    std::vector<int> neighbors_index;//存储邻居索引号
    std::vector<double> neighbor_cost;//存储选择每个邻居作为父节点的消耗

    int bestparent_index = 0;
    double bestparent_dis = 10000;
    int neighbor_num = 0;//邻居个数
    double neighbor2new_dis = 0.0;
    int length = (int)trees.size();
    VectorXd temp;//没啥用，存储返回的扩展点
    for(int i = 0; i < length; i++ )
    {
        temp_differ = vector_new - trees[i].vector;
        neighbor2new_dis = temp_differ.norm();
        if(neighbor2new_dis < radius)
        {
            //首先尝试连接第i个节点和vector_new，若在中途碰撞，则舍弃第i个节点
            connect_neighber_result = connectCSpaceRRTStar(trees, vector_new, temp, i);
            if(connect_neighber_result == connection_success)
            {
                neighbors_index.push_back(i) ;//保存邻居编号
                neighbor_cost.push_back(trees[i].dis + neighbor2new_dis) ;//保存所选父节点的cost
                //是否优化距离了？
                if(neighbor_cost[neighbor_num] < bestparent_dis)
                {
                    bestparent_dis = neighbor_cost[neighbor_num];
                    bestparent_index = i;//保存最优父节点索引号
                }
                neighbor_num++;
            }
        }
    }
    //除去最优父节点，返回剩余的周围节点
//    int k = 0;
    for(int i=0; i < (int)neighbors_index.size(); i++)
    {
        if(neighbors_index[i] != bestparent_index)
        {
            neighbors_others.push_back(neighbors_index[i]);
        }
    }

    return bestparent_index;
}

/** rewire之后，需要修改子节点的cost
 *  输入：
 *       trees：                     tree_start 或 tree_goal
 *       nodeIndex：                 尝试进行子节点cost更新的节点
 *       rewiredNeighborIndex：      rewire后得到优化的Neighbor的Index
 *  输出：
 *       直接在trees中更新各子节点的cost
 *
 * */
void RRTStarPlanner::modifyChildCost(std::vector<Path_Node> trees, int nodeIndex, const int rewiredNeighborIndex)
{
    //退出递归的条件， 如果该节点child为空,则更新该节点cost
    if(trees[nodeIndex].child.empty())
    {
        VectorXd differTemp;
        const int nodeIndexTemp = nodeIndex;//临时保存index
        int parentIndex;//父节点index
        double cost = 0.0;
        while(nodeIndex != rewiredNeighborIndex)
        {
            parentIndex = trees[nodeIndex].parent;
            differTemp = trees[parentIndex].vector - trees[nodeIndex].vector;
            cost = cost + differTemp.norm();
            nodeIndex = parentIndex;//更新Index,往rewireNeighborIndex逼近
        }
        cost = cost + trees[rewiredNeighborIndex].dis;
        trees[nodeIndexTemp].dis = cost;//更新该节点cost

    } else
    {
        int len = (int)trees[nodeIndex].child.size();//子节点个数
        for(int i = 0; i < len; i++)
        {
            //           std::cout << trees[nodeIndex].child[i] << std::endl;//
            modifyChildCost(trees, trees[nodeIndex].child[i], rewiredNeighborIndex);//递归
        }
    }
}