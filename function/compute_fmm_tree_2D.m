function tree = compute_fmm_tree_2D(unique_coords, unicoord_indics, method, max_points, min_points, eta, cheb_num, A, A_x,concentration_threshold,ar_thresh,use_parallel)
    % COMPUTE_FMM_TREE 使用 FMM 计算树结构
    % 输入参数:
    %   unique_coords - 二维点集 [N×2]，每行是一个点的 x,y 坐标
    %   unicoord_indics - 点集的索引集合 [N×1]，与 unique_coords 对应
    %   method - 旧划分方法的四叉树划分方式（1 或 2）
    %   max_points - 每个盒子最大点数
    %   min_points - 每个盒子最小点数
    %   eta - 分离系数，控制不可接受区域范围
    %   cheb_num - 一维方向切比雪夫节点数量
    %   A - 核函数句柄 A(x, y)，输入两个二维点 x, y，返回标量
    %   A_x - 句柄函数，接受四个向量(x1,x2,y1,y2)，返回矩阵，用于向量化计算
    %   use_parallel - 是否启用并行计算（0 或 1）
    % 输出:
    %   tree - 最终的 FMM 树结构（对应原代码中的 tree_5）

    % 构造映射函数句柄 S(a,b)
    % S = construct_Sab_chebfun_2D(cheb_num); % a 是标准区间上的切比雪夫插值点，b 是标准区间上的函数值点
    
    % 向量化转移函数
    S_x = construct_Sab_chebfun_x_2D(cheb_num); % 向量化转移函数


    % ========= 旧划分方法 =================== %
    % % 构建四叉树
    % tree_1 = quadtree_2D_decomposition(unique_coords, unicoord_indics, ...
    %     max_points, min_points, method,concentration_threshold);
    % 
    % % 可接受性判断与树结构重构
    % tree_2 = analyze_level_admissibility_2D(tree_1, eta, use_parallel);

    % ========= 新划分方法（混合树） ========== %

    % 构建混合四叉树(直接重构)
    [~, level_nodes_struct] = mixquadtree_2D_decomposition(unique_coords, unicoord_indics, ...
        max_points, min_points,concentration_threshold,ar_thresh);

    % 可接受性判断
    tree_2 = analyze_level_admissibility_mixtree_2D(level_nodes_struct, eta, use_parallel);
    
    % ====================================== %

    % 使用 chebpts 生成一维切比雪夫节点（第一类切比雪夫点，不包含端点）
    cheb_p = chebpts(cheb_num, 1);

    % 张量积生成二维网格节点
    [X_c, Y_c] = meshgrid(cheb_p, cheb_p);

    % 将节点转换为列向量形式，y 优先排序
    cheb_nodes = sortrows([X_c(:), Y_c(:)], [2, 1]);

    % 将标准区间切比雪夫节点映射到物理单元的每一层的每个盒子上
    tree_3 = map_to_box_elements_2D(cheb_nodes, tree_2, use_parallel);

    % 将物理单元的每个盒子中的核函数点集映射到上一层的标准区间上
    tree_4 = map_points_2D(tree_3, use_parallel);

   % 同时计算整层和各个盒子的转移矩阵和系数矩阵
   % 向量化计算
   tree_5 = compute_M2M_L2L_A_x_2D(S_x, cheb_nodes, A_x, tree_4, length(unique_coords), unique_coords, use_parallel);
    
   % 非向量化计算
   % tree_5 = compute_M2M_L2L_A_2D(S, cheb_nodes, A, tree_4, length(unique_coords), unique_coords, use_parallel);

    % 返回最终树结构
    tree = tree_5;
end