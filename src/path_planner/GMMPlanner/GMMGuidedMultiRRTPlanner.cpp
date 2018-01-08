//
// Created by zhaoxin on 17-12-28.
//

#include "GMMGuidedMultiRRTPlanner.hpp"

void MultiRRTTree::merge(MultiRRTTree &mergedTree, int newParentIndex, int oldIndex)
{
    // 先将树编号添加进来
    for(int i = 0; i < mergedTree.index.size(); i++)
        index.push_back(mergedTree.index[i]);

    mergedTree.isMerged = true;

    // 是否已经添加的标志位
    bool *isAddedFlag;
    isAddedFlag = new bool[mergedTree.pathNode_tree.size()]();

    // 将连接点的子节点都添加进来
    int newIndex = (int)pathNode_tree.size();
    addSubTree(mergedTree, newParentIndex, oldIndex, isAddedFlag);

    isAddedFlag[oldIndex] = true;

    int oldParentIndex;
    // 循环添加
    do{
        oldParentIndex = mergedTree.pathNode_tree[oldIndex].parent;
        oldIndex = oldParentIndex;
        newIndex = addSubTree(mergedTree, newIndex, oldParentIndex, isAddedFlag);
    }while(oldParentIndex != 0);

    delete [] isAddedFlag;
}

/**
 * 将该节点的其他子节点及其分支都添加进来
 * @param mergedTree     ： 要合并的树
 * @param newParentIndex ： 此节点在新树中的父节点编号
 * @param oldParent      ： 在待合并的树中的
 * @param addedNode
 */
int MultiRRTTree::addSubTree(MultiRRTTree mergedTree, int newParentIndex, int currentOldIndex, bool *isAddedFlag)
{
    // 计算距离
    VectorXd vector_temp = pathNode_tree[newParentIndex].vector - mergedTree.pathNode_tree[currentOldIndex].vector;
    // 将当前节点添加到新树上去
    pathNode_tree.push_back(Path_Node(mergedTree.pathNode_tree[currentOldIndex].vector, newParentIndex, vector_temp.norm()));

    // 在新树上的编号
    int newIndex = (int)pathNode_tree.size() - 1;

    // 将当前新的编号添加到父节点的子树中去
    pathNode_tree[newParentIndex].child.push_back(newIndex);

    // 标记此节点已经添加了
    isAddedFlag[currentOldIndex] = true;

    for(int i = 0; i < mergedTree.pathNode_tree[currentOldIndex].child.size(); i++)
    {
        // 如果这个子节点并没有被添加， 则迭代添加子节点
        if(!isAddedFlag[mergedTree.pathNode_tree[currentOldIndex].child[i]])
        {
            // 迭代添加所有节点
            addSubTree(mergedTree, newIndex, mergedTree.pathNode_tree[currentOldIndex].child[i], isAddedFlag);
        }
    }

    // 返回进入此程序时的父节点编号值
    return newIndex;
}

/**
 * 初始化规划器
 * @param str  ： 起始位置
 * @param goal ： 目标位置
 */
void GMMGuidedMultiRRTPlanner::Init(VectorXd str, VectorXd goal)
{
    _str = str;
    _goal = goal;

    step = 0.01; // 碰撞检测的步长

    // 关节采样范围
    randomRangeLow.resize(6);
    randomRangeUp.resize(6);

    randomRangeLow << -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI;
    randomRangeUp << M_PI, M_PI, M_PI, M_PI, M_PI, M_PI;

    trees.resize(GMM.CSpaceGmm->nstates);

    // 在每个中心处初始化
    VectorXd center_vector_temp;
    center_vector_temp.resize(6);
    float *center_arr_temp = new float[6]();
    for(int i = 0; i < GMM.CSpaceGmm->nstates; i++)
    {
        GMM.CSpaceGmm->GetMean(i, center_arr_temp);
        arr2vector(center_arr_temp, center_vector_temp);
        trees[i].pathNode_tree.push_back(Path_Node(center_vector_temp, 0, 0));
        trees[i].isMerged = false;
        trees[i].index.push_back(i);
    }

    E.resize(GMM.CSpaceGmm->nstates, GMM.CSpaceGmm->nstates);
    disMat.resize(GMM.CSpaceGmm->nstates, GMM.CSpaceGmm->nstates);

    vexNum = GMM.CSpaceGmm->nstates;

    threshold_sampling = 0.2;
    threshold_greedy = 0.8;
}

PathSearch_result GMMGuidedMultiRRTPlanner::plan()
{

    strGaussionIndex = findNearestGaussian(_str);
    goalGaussionIndex = findNearestGaussian(_goal);

    PathSearch_result pathSearchResult;
//
//    std::vector<Path_Node> strTree, goalTree;

    std::vector<VectorXd> startPath = planEndPath(getCenter(strGaussionIndex), _str);
    std::vector<VectorXd> goalPath = planEndPath(getCenter(goalGaussionIndex), _goal);
    path.clear();

    std::cout << "start path size:" << startPath.size() << "goal path: " << goalPath.size() << std::endl;

    //第一步，先规划str到第strIndex中心点的路径，此路径采用标准的单向RRT算法
//    pathSearchResult = singleRRTSearch(_str, strGaussionIndex, strTree);
//    if(pathSearchResult == pathsearch_fail)
//        return pathsearch_fail;
//
//    std::vector<VectorXd> plannedPathStr;
//    //将第一段路径添加到已规划好的路径中
//    std::vector<VectorXd>::iterator it;
//    int insertIndex = (int)strTree.size() - 1;
//    while(insertIndex != 0)
//    {
//        it = plannedPathStr.begin();
//        plannedPathStr.insert(it, strTree[insertIndex].vector);
//        insertIndex = strTree[insertIndex].parent;
//    }
//    it = plannedPathStr.begin();
//    plannedPathStr.insert(it, strTree[0].vector);
//    path.clear();
//    path = plannedPathStr;
//    simplifyPath();
//    plannedPathStr = path;


//    trees[strGaussionIndex].pathNode_tree = tree_start;

    // 第二步，从目标点规划到最近
//    pathSearchResult = singleRRTSearch(_goal, goalGaussionIndex, goalTree);
//    if(pathSearchResult == pathsearch_fail)
//        return pathsearch_fail;
//    trees[goalGaussionIndex].pathNode_tree = goalTree;

    // 第三步，规划strIndex到goalIndex中心的路径，此部分采用双向RRT算法，并且在利用PRM规划出来的最优路径附近采样。
    preV = findPathDigkstra(strGaussionIndex);    // 找到最优的中心路径

    int addvexIndex = goalGaussionIndex;

    sampleIndexVector.push_back(addvexIndex);

    std::cout << "start gaussion index: " << strGaussionIndex << std::endl;
    std::cout << "gual gaussion index: " << goalGaussionIndex << std::endl;


    while(addvexIndex != strGaussionIndex)
    {
        addvexIndex = preV[addvexIndex];
        sampleIndexVector.push_back(addvexIndex);

    }

    std::cout << "GMM states: ";
    for(int i = 0; i < sampleIndexVector.size(); i++)
        std::cout << sampleIndexVector[i] << ", ";
    std::cout << std::endl;

    PathSearch_result resultE = multiRRTExtend();


    std::vector<VectorXd>::iterator it;
    if(resultE == pathsearch_finish)
    {
        for(int i = 0; i < startPath.size(); i++)
        {
            it = path.begin();
            path.insert(it, startPath[i]);
        }

        for(int i = 0; i < goalPath.size(); i++)
        {
            path.push_back(goalPath[i]);
        }
        simplifyPath();
        return pathsearch_finish;
    }
    else{
        return pathsearch_fail;
    }

}

//std::vector<int> sampleIndexVector
PathSearch_result GMMGuidedMultiRRTPlanner::multiRRTExtend()
{
    double epsilon;

    ExtendTree_result extendResult = extend_start_fail, subExtandResult = extend_start_fail;
    int sampleIndex1, sampleIndex2, sampleIndexTreeIndex, sampleTreeIndex_1, sampleTreeIndex_2;

    int strTreeIndex, goalTreeIndex;
    int extendNum = 0;
    while(extendResult == extend_start_fail && extendNum < 2000000)
    {
        epsilon = getrandom(0., 1.); // 生成随机数，根据随机数大小进行选择
        if(epsilon < threshold_sampling)
        {
            std::cout << "random sampling" << std::endl;
            // 按照一定概率在空间中随机采样
            subExtandResult = extendBiRRT(strGaussionIndex, goalGaussionIndex, false);
            if(subExtandResult == extend_finish)
                extendResult = extend_finish;
        } else
        {
            goalTreeIndex = findGaussianIndexInTrees(goalGaussionIndex);
            // 否则在GMM内部采样
            if(epsilon < threshold_greedy)
            {
                std::cout << "greedy sampling" << std::endl;
                // 采用greedy策略，直接找Tree_i中离goal最近的高斯G_i与下一个高斯G_{i+1}进行RRT操作
                sampleIndex1 = (int)random();
                sampleIndex1 = sampleIndex1 % (int)sampleIndexVector.size();
                sampleTreeIndex_1 = findGaussianIndexInTrees(sampleIndexVector[sampleIndex1]);
                if(sampleTreeIndex_1 == goalTreeIndex)
                    continue;

                sampleIndex2 = findGreedyIndex(sampleIndex1);
                sampleTreeIndex_2 = findGaussianIndexInTrees(sampleIndex2);

            }
            else
            {
                std::cout << "GMM sampling" << std::endl;
                // 否则，在选择不在Tree_i中的其他G_k，进行G_i与G_i的双向RRT操作
                sampleIndex1 = (int)random();
                sampleIndex1 = sampleIndex1 % GMM.CSpaceGmm->nstates;
                sampleIndex2 = (int)random();
                sampleIndex2 = sampleIndex2 % GMM.CSpaceGmm->nstates;
                sampleTreeIndex_1 = findGaussianIndexInTrees(sampleIndex1);
                sampleTreeIndex_2 = findGaussianIndexInTrees(sampleIndex2);
            }
            if(sampleTreeIndex_1 == sampleTreeIndex_2)
                continue;
            extendBiRRT(sampleTreeIndex_1, sampleTreeIndex_2, true);
        }

        strTreeIndex = findGaussianIndexInTrees(strGaussionIndex);
        goalTreeIndex = findGaussianIndexInTrees(goalGaussionIndex);
        if(strTreeIndex == goalTreeIndex)
            extendResult = extend_finish;
    }

    int samplenum = 0;
    for(int i = 0; i < trees.size(); i++)
    {
        if(!trees[i].isMerged)
        {
            samplenum += (int)trees[i].pathNode_tree.size();
        }
    }
    std::cout << "total sample number is: " << samplenum << std::endl;

    return extendResult == extend_finish ? pathsearch_finish : pathsearch_fail;
}

/**
 * 采用greedy策略找到下一个最优的tree编号
 * @param currentGaussionIndex
 * @return
 */
int GMMGuidedMultiRRTPlanner::findGreedyIndex(int currentGaussionIndex)
{
    // 先在最优路径里面找到currentIndex所对应的最优高斯序列里的编号
    int sampleVectIndex = 0;
    for(int i = 0; i < sampleIndexVector.size(); i++)
    {
        if(sampleIndexVector[i] == currentGaussionIndex)
            sampleVectIndex = i;
    }

    int greedyGaussionIndex;

    // 判断currentIndex是否和goalIndex在同一个tree中，如果在一个树中则退出，并返回2018表示结束
    int goalTreeIndex = findGaussianIndexInTrees(goalGaussionIndex);
    int currentTreeIndex = findGaussianIndexInTrees(currentGaussionIndex);
    if(goalTreeIndex == currentTreeIndex)
        return 2018;

    // TODO: 这里存在问题
    // 否则找到下一个并与currentIndex不在同一个tree中的index

    int greedyTreeIndex;
    for(int i = 0; i < sampleIndexVector.size(); i++)
    {
        greedyGaussionIndex = sampleIndexVector[i];
        greedyTreeIndex = findGaussianIndexInTrees(greedyGaussionIndex);
        if(greedyTreeIndex != currentTreeIndex)
            break;
    }

//    do{
//        greedyGaussionIndex = sampleIndexVector[++sampleVectIndex];
//        greedyTreeIndex = findGaussianIndexInTrees(greedyGaussionIndex);
//    }while(greedyTreeIndex == currentTreeIndex);

    return greedyTreeIndex;
}

/**
 * 找到第gaussianIndex 个高斯所在的 tree 的 index
 * @param gaussianIndex : 要找到的高斯索引值
 * @return： 所在的tree的索引值
 */
int GMMGuidedMultiRRTPlanner::findGaussianIndexInTrees(int gaussianIndex)
{
    for(int i = 0; i < trees.size(); i++)
    {
        if(!trees[i].isMerged)
        {
            if(std::find(trees[i].index.begin(), trees[i].index.end(), gaussianIndex) != trees[i].index.end())
            {
                return i;
            }
        }
    }
}

/**
 * 在第i个高斯和第j个高斯之间采用双向RRT进行规划
 * @param i
 * @param k
 * @param isInGmm : 标志位，如果为真，则在GMM内部生成随机点，否则随机生成
 * @return
 */
ExtendTree_result GMMGuidedMultiRRTPlanner::extendBiRRT(int sample1, int sample2, bool isInGmm)
{
    std::cout << "extend trees:" << sample1 << " " << sample2 << std::endl;
    int randi = (int)random();
    randi = randi % 2;

    int sampleIndex = (randi == 0) ? sample1 : sample2;

    VectorXd strVector = getCenter(sample1);
    VectorXd goalVector = getCenter(sample2);

    // 创建一个保存连接树的向量
    VectorXd vector_out(6);

    Connection_result connection_result;

    // 随机生成一个无碰撞点
    float randomArr[6];
    bool collision_flag;
    VectorXd vector_random(6);
    do{
        // 这里根据选项是否在GMM内采样进行处理
        if(isInGmm)
        {
            // 如果在GMM内部规划，那么随机点利用GMM模型生成
            std::cout << "sample state: " << sampleIndex << std::endl;
            gaussian_draw(&(GMM.CSpaceGmm->c_gmm->gauss[sampleIndex]), randomArr);
            GMM.CSpaceGmm->Pdf(randomArr, sampleIndex);
            arr2vector(randomArr, vector_random);       // 转换到vector中
        } else{
            // 否则随机生成
            for(int i = 0; i < 6; i++)
                vector_random[i] = getrandom(randomRangeLow[i], randomRangeUp[i]);
        }
//        std::cout << "random vector" <<vector_random.transpose() << std::endl;
        collision_flag = collisionCheckRight(vector_random);    // 判断此点是否存在碰撞， 如果碰撞，那么继续生成
    }while(collision_flag);

    VectorXd vector_temp_str(6), vector_temp_goal(6);
    // 尝试进行连接
    connection_result = connectCSpaceRRT(trees[sample1].pathNode_tree, vector_random, vector_temp_str);
    if(connection_result == connection_fail)
    {
//        ROS_INFO("extend start tree failed");
        return extend_start_fail;
    } // 如果不能连接，则返回失败
    else if(connection_result == connection_success)
        vector_temp_str = vector_random;    // 如果连接成功，则将最终连接的点设为这个随机点

//    ROS_INFO("extend start tree success");
    connection_result = connectCSpaceRRT(trees[sample2].pathNode_tree, vector_temp_str, vector_temp_goal);

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

        if(trees[sample1].index[0] == strGaussionIndex)
        {
            if(trees[sample2].index[0] == goalGaussionIndex){
                pathFinish(sample1, sample2); }

            trees[sample1].merge(trees[sample2], (int)trees[sample1].pathNode_tree.size() - 1, (int)trees[sample2].pathNode_tree.size() - 1);
        } else if(trees[sample1].index[0] == goalGaussionIndex)
        {
            if(trees[sample2].index[0] == strGaussionIndex)
            {
                pathFinish(sample2, sample1);
                trees[sample2].merge(trees[sample1], (int)trees[sample2].pathNode_tree.size() - 1, (int)trees[sample1].pathNode_tree.size() - 1);// 合并tree_i和tree_k
            }
            else{
                trees[sample1].merge(trees[sample2], (int)trees[sample1].pathNode_tree.size() - 1, (int)trees[sample2].pathNode_tree.size() - 1);
            }
        }
        else{
            // 这里面包含三种情况
            trees[sample2].merge(trees[sample1], (int)trees[sample2].pathNode_tree.size() - 1, (int)trees[sample1].pathNode_tree.size() - 1);// 合并tree_i和tree_k
        }

        // TODO: 输出当前树结构数据，方便查看，后面可以取消
        for(int i = 0; i < trees.size(); i++)
        {
            if(!trees[i].isMerged)
            {
                std::cout << "第 " << i << "个树，包含节点：";
                for(int j = 0; j < trees[i].index.size(); j++)
                {
                    std::cout << trees[i].index[j] << ", ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << std::endl << std::endl;
        vector_temp_goal = vector_temp_str;
    }
    else
    {
//        ROS_INFO("extend success");
        extend_result = extend_success;
    } // 如果连接上了一部分，则返回部分成功

    return extend_result;
}

void GMMGuidedMultiRRTPlanner::pathFinish(int str, int goal)
{

    // 先清除path
    path.clear();

    // 将规划的树重新排序到输出的变量rrt_path中去，并返回规划成功
    std::vector<VectorXd>::iterator it;

    int i = (int)trees[str].pathNode_tree.size() - 1;

    // 将tree_start逆向排列
    while(i != 0)
    {
        it = path.begin();
        path.insert(it, trees[str].pathNode_tree[i].vector);
        i = trees[str].pathNode_tree[i].parent;
    }
    it = path.begin();
    path.insert(it, trees[str].pathNode_tree[0].vector);

    // 逐个添加tree_goal上的节点
    i = (int)trees[goal].pathNode_tree.size() - 1;

    while(i !=0 )
    {
        path.push_back(trees[goal].pathNode_tree[i].vector);
        i = trees[goal].pathNode_tree[i].parent;
    }
    path.push_back(trees[goal].pathNode_tree[0].vector);
}
