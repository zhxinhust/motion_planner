//
// Created by zhaoxin on 17-12-28.
//

#ifndef MOTION_PLANNERS_GMMGUIDEDMULTIRRTPLANNER_HPP
#define MOTION_PLANNERS_GMMGUIDEDMULTIRRTPLANNER_HPP

#include "GMMGuidedPlanner.hpp"

class GMMGuidedMultiRRTPlanner:public GMMGuidedPlanner
{
public:
    GMMGuidedMultiRRTPlanner(robot_state::RobotState &a, Robot_moveit &robMoveit):GMMGuidedPlanner(a, robMoveit)
    {

    }

    void Init(VectorXd str, VectorXd goal);

    PathSearch_result plan();

    PathSearch_result multiRRTExtend();

    ExtendTree_result extendBiRRT(int i, int k, bool isInGmm);

    int findGreedyIndex(int currentIndex);

    int findGaussianIndexInTrees(int gaussianIndex);

    std::vector<MultiRRTTree> trees;

protected:
    // 根据生成的随机数ε来选择采样点策略，如果ε < threshold_sampling 则在运动空间中随机生成采样点，
    // 如果threshold_sampling < ε < threshold_greedy 则采用greedy策略，搜索最优路径策略中，不在Tree_i中的下一个最优路径点
    // 如果 threshold_greedy < ε 则 随机找一个不在Tree_i中的Tree
    double threshold_greedy;
    double threshold_sampling;

    int strGaussionIndex;   // 离起始点最近的
    int goalGaussionIndex;

    VectorXi preV;
    std::vector<int> sampleIndexVector;
};

class MultiRRTTree
{
public:
    MultiRRTTree() {}

    std::vector<int> index; // 包含的中心数，如果两个RRT树相连了，则进行合并处理，这样一个树就包含有多个中心了
    std::vector<Path_Node> pathNode_tree;    // 内部包含有树节点，每个节点上有距离、父节点和位置信息

    bool isMerged;  // 此变量标志是否已经被合并，如果被合并了，则跳过这个结构的搜索

    void merge(MultiRRTTree tree, int parentIndex, int childIndex);

    int addSubTree(MultiRRTTree mergedTree, int newParentIndex, int currentIndex, bool *isAddedFlag);
};

#endif //MOTION_PLANNERS_GMMGUIDEDMULTIRRTPLANNER_HPP
