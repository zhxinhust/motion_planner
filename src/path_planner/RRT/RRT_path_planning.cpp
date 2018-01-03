//
// Created by zhxin on 17-8-3.
//


//#include "RRT_path_planning.h"
#include "RRTPlanner.h"

/**
 * 单臂RRT规划算法
 * @param str
 * @param goal
 * @param robot_moveit
 * @param rrt_path
 * @param pub
 * @return
 */
PathSearch_result rrt_planning(VectorXd str, VectorXd goal, Robot_moveit *robot_moveit,
                               std::vector<VectorXd> &rrt_path, ros::Publisher pub)
{
    std::vector<Path_Node> tree_start, tree_end;    // 建立处于起点和目标点处的tree结构

    Path_Node pathnode_temp(6);    // 采样用的path node

    // 先创建start处的树结构
    pathnode_temp.vector = str;
    tree_start.push_back(pathnode_temp);

    // 创建goal 处的树结构
    pathnode_temp.vector = goal;
    tree_end.push_back(pathnode_temp);

    // 创建一个保存连接树的向量
    VectorXd vector_out;
    vector_out.resize(6);

    // 获取创建碰撞检测所需要的变量
    robot_state::RobotState &check_collision_state = robot_moveit->planning_scene->getCurrentStateNonConst();

    // 先尝试直接连接start与goal，如果中间无碰撞，则说明这两个点可以直接相连接，那么就不需要其他规划了。
    // Collision_result collision_result = cspace_connection(str, goal, robot_moveit, vector_out, check_collision_state, pub);

    Collision_result collision_result = cspace_connection(tree_start, 0, goal, robot_moveit,
                                                          vector_out, check_collision_state, pub);
    // 如果中间无碰撞，则可以直接退出了
    if(collision_result == collision_free)
    {
        rrt_path.push_back(str);
        rrt_path.push_back(goal);

        return pathsearch_finish;
    }

    int extendnum = 0;
    ExtendTree_result extend_result = extend_start_fail;

    // 一直扩展树，直到到达预定次数，或者规划成功
    while(extendnum < 10000 && extend_result != extend_finish)
    {
        // 扩展树
        extend_result = cspace_extendTree(tree_start, tree_end, robot_moveit, check_collision_state, pub);
        extendnum++;
    }

    // 如果到达了规定的次数还没有规划成功，则返回规划失败
    if(extendnum == 10000)
        return pathsearch_fail;

    // 将规划的树重新排序到输出的变量rrt_path中去，并返回规划成功
    std::vector<VectorXd>::iterator it;

    std::cout << "start tree size:" << tree_start.size() << std::endl;
    std::cout << "end tree size:" << tree_end.size() << std::endl;


    int i = (int)tree_start.size() - 1;

    while(i != 0)
    {
        it = rrt_path.begin();
        rrt_path.insert(it, tree_start[i].vector);
        i = tree_start[i].parent;
    }
    it = rrt_path.begin();
    rrt_path.insert(it, tree_start[0].vector);

    i = (int)tree_end.size() - 1;
    while(i !=0 )
    {
        rrt_path.push_back(tree_end[i].vector);
        i = tree_end[i].parent;
    }
    rrt_path.push_back(tree_end[0].vector);

    return pathsearch_finish;
}

/**
 * C 空间中，连接str 与 goal进行碰撞检测，如果
 * @param str
 * @param goal
 * @param robot_moveit
 * @param vector_out
 * @param check_collision_state
 * @return 如果两者之间不存在碰撞，则返回0,如果第一步就碰撞，则返回1, 如果是路径中间才发生碰撞，则返回2,并且将最后无碰撞的点保存在vector_out中
 */

Collision_result cspace_connection(std::vector<Path_Node> &trees, int nearestindex, VectorXd goal,
                                   Robot_moveit *robot_moveit, VectorXd &vector_out,
                                   robot_state::RobotState &check_collision_state, ros::Publisher pub)
{

    VectorXd str = trees[nearestindex].vector;

    // 输出信息
    std::cout << "try to connect str:";
    for(int i = 0; i < 6; i++)
        std::cout<< str[i] << ' ';
    std::cout << "goal: ";

    for(int i = 0; i < 6; i++)
        std::cout << goal[i] << ' ';
    std::cout << std::endl;

    const double step = 0.05;   // 碰撞检测的步长

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
        add_tree_node(trees, nearestindex, goal);
        return collision_free;
    }

    // 以最小步长进行尝试连接，直到连接到了目标点，或者发生碰撞
    while((!iscollision) && len < vector_temp.norm())
    {
        len = len + step;
        vector_out = str + len * vectorUnitDirection;
        iscollision = collision_detection_right(vector_out, robot_moveit, check_collision_state, pub);

        // 如果不碰撞，则添加到树上去
        if(!iscollision)
        {
            parent_index = (addnum == 0) ? nearestindex : (int)trees.size() - 1;
            add_tree_node(trees, parent_index, vector_out);

            addnum++;
        }
    }

    // 如果两节点不发生碰撞，则返回0
    if (vector_temp.norm() - len < step)
    {
        ROS_INFO("connection is free");
        return collision_free;
    }
    else if(len <= step)
    {
        ROS_INFO("connection is invalid");
        return collision_invalid;
    }
    else
    {
        ROS_INFO("A mid point is found");
        // 最后一次保存的是碰撞的点，这里减去增量，得到最后无碰撞的点
        vector_out = vector_out - step * vectorUnitDirection;
        return collision_mid;
    }
}

/**
 * 给树添加节点
 * @param trees  ：树
 * @param parent ：父节点下标
 * @param v      ：添加的节点向量
 */
void add_tree_node(std::vector<Path_Node> &trees, int parent, VectorXd v)
{
    Path_Node path_node_new;

    path_node_new.vector = v;
    path_node_new.parent = parent;

    VectorXd v_temp = v - trees[parent].vector;
    path_node_new.dis = v_temp.norm() + trees[parent].dis;

    trees.push_back(path_node_new);
}

ExtendTree_result cspace_extendTree(std::vector<Path_Node> &str_tree, std::vector<Path_Node> &goal_tree,
                                    Robot_moveit *robot_moveit, robot_state::RobotState &check_collision_state,
                                    ros::Publisher pub)
{
    VectorXd vector_random, vector_temp_str, vector_temp_goal, dis_temp;
    Collision_result collision_result;
    Path_Node path_node_new(6);

    bool collision_flag = true;

    // 各关节的采样范围
    double random_low[] = {-M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI};
    double random_up[] = {M_PI, M_PI, M_PI, M_PI, M_PI, M_PI};

    vector_random.resize(6);
    vector_temp_str.resize(6);
    vector_temp_goal.resize(6);

    // 随机生成一个无碰撞的点
    while(collision_flag)
    {
        for(int i = 0; i < 6; i++)
            vector_random[i] = getrandom(random_low[i], random_up[i]);

        collision_flag = collision_detection_right(vector_random, robot_moveit, check_collision_state, pub);
    }

    int nearestindex = cspace_findNearestNode(vector_random, str_tree); // 在start树上找到距离最近的点

    // 尝试进行连接
  //  collision_result = cspace_connection(str_tree[nearestindex].vector,
  //                                       vector_random, robot_moveit, vector_temp_str, check_collision_state, pub);

    collision_result = cspace_connection(str_tree, nearestindex, vector_random,
                                         robot_moveit, vector_temp_str,check_collision_state, pub);

    if(collision_result == collision_invalid)
    {
        ROS_INFO("extend start tree failed");
        return extend_start_fail;
    } // 如果不能连接，则返回失败
    else if(collision_result == collision_free)
        vector_temp_str = vector_random;    // 如果连接成功，则将最终连接的点设为这个随机点

    ROS_INFO("extend start tree success");

    // 找到与刚添加到start树上的节点距离最近的点
    nearestindex = cspace_findNearestNode(vector_temp_str, goal_tree);
    // 尝试进行连接
   // collision_result = cspace_connection(goal_tree[nearestindex].vector,
   //                                      vector_temp_str, robot_moveit, vector_temp_goal, check_collision_state, pub);

    collision_result = cspace_connection(goal_tree, nearestindex, vector_temp_str,
                                         robot_moveit, vector_temp_str,check_collision_state, pub);

    ExtendTree_result extend_result;

    if(collision_result == collision_invalid)
    {
        ROS_INFO("extend end fail");
        return extend_end_fail;
    }     // 如果连接不成功，则返回失败
    else if (collision_result == collision_free)
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


/**
 * 在C空间中找到树上最近的点
 * @param p
 * @param trees
 * @return 返回树上节点的编号
 */
int cspace_findNearestNode(VectorXd p, std::vector<Path_Node> trees)
{
    double dis = 10000;
    int index = 0;
    VectorXd temp;

    for(unsigned int i = 0; i < trees.size(); i++)
    {
        temp = trees[i].vector - p;
        if(temp.norm() < dis)
        {
            index = i;
            dis = temp.norm();
        }
    }

    printf_pos("nearest point is: ", trees[index].vector);
    //std::cout << "nearest is " << index << ": " << trees[index].vector << std::endl;
    return index;
}

/**
 * 检测右臂的碰撞结果
 * @param p             ： 右臂关节位置
 * @param robot_moveit  ： moveit 类
 * @param current_state
 * @return              ： 如果碰撞，则返回1，否则返回0
 */
bool collision_detection_right(VectorXd p, Robot_moveit *robot_moveit, robot_state::RobotState & current_state, ros::Publisher pub)
{
    // 将要检测的角度位置更新到模型中去
    std::vector<double> joint_pos(6, 0);

    // 建立jointstate变量，用来发布关节信息
    JointState joint = InitJointState(JntArray(16));

    for(unsigned i = 0; i < 6; i++)
    {
        joint_pos[i] = p[i];
        joint.position[i + 8] = p[i];
    }

    // 发布关节位置
    joint.header.stamp = ros::Time::now();
    pub.publish(joint);

    current_state.setJointGroupPositions(robot_moveit->joint_model_group_right, joint_pos);

    // 清除原来计算结果
    robot_moveit->collision_result_right.clear();

    // 碰撞检测
    robot_moveit->planning_scene->checkCollision(*robot_moveit->collision_request_right,
                                                robot_moveit->collision_result_right);

    // 输出现在检测的碰撞位置及结果
    //printf_pos("current collision detect pos:", p);
    std::cout << "current collision detect pos:" ;
    for(int i = 0; i < 6; i++)
    {
        std::cout<< p[i] << ' ';
    }
    std::cout<<(robot_moveit->collision_result_right.collision ? "is collision" : "is free") << std::endl;

    // 返回碰撞检测结果
    return robot_moveit->collision_result_right.collision;
}



/**
 * 输出角度及变量，调试用
 * @param str
 * @param p
 */
void printf_pos(const char *str, VectorXd p)
{
    std::cout << str;
    for(long i = 0; i < p.size(); i++)
        std::cout << p[i] << ' ';
    std::cout << std::endl;
}