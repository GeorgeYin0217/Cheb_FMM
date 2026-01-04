function mapped_struct_tree = assemble_layer_M2M_L2L_A_2D(mapped_struct_tree, use_parallel)
    % 输入：
    % mapped_struct_tree: 已包含每盒子M2M、L2L、A矩阵的结构数组
    % use_parallel: 0表示串行，1表示并行计算
    %
    % 输出：
    % mapped_struct_tree: 更新后的结构数组，新增字段：
    %   - M2M_layer: 每层所有盒子的M2M矩阵组装成的块对角矩阵
    %   - L2L_layer: 每层所有盒子的L2L矩阵组装成的块对角矩阵
    %   - A_layer: 每层所有盒子的A矩阵组装成的块对角矩阵
    %
    % 索引规则：
    % 假设每层有n_boxes个子盒，M2M为n_cheb×n_map，L2L为n_map×n_cheb，
    % A为n_cheb_local×n_neighbor（n_neighbor = n_neighbors_per_box × n_cheb_local）。
    %
    % 1. M2M_layer：
    %    - 行索引：总行数为n_boxes×n_cheb，按盒子顺序分块排列
    %      - 盒子k的行：(k-1)×n_cheb + 1 到 k×n_cheb
    %      - 示例：n_boxes=3, n_cheb=16，总行数48
    %        盒子1：1–16, 盒子2：17–32, 盒子3：33–48
    %    - 列索引：总列数为n_boxes×n_map，按盒子顺序分块排列
    %      - 盒子k的列：(k-1)×n_map + 1 到 k×n_map
    %      - 示例：n_boxes=3, n_map=20，总列数60
    %        盒子1：1–20, 盒子2：21–40, 盒子3：41–60
    %    - 结构：块对角矩阵 [M2M_1, 0, 0; 0, M2M_2, 0; 0, 0, M2M_3]
    %
    % 2. L2L_layer：
    %    - 行索引：总行数为n_boxes×n_map，按盒子顺序分块排列
    %      - 盒子k的行：(k-1)×n_map + 1 到 k×n_map
    %      - 示例：n_boxes=3, n_map=20，总行数60
    %        盒子1：1–20, 盒子2：21–40, 盒子3：41–60
    %    - 列索引：总列数为n_boxes×n_cheb，按盒子顺序分块排列
    %      - 盒子k的列：(k-1)×n_cheb + 1 到 k×n_cheb
    %      - 示例：n_boxes=3, n_cheb=16，总列数48
    %        盒子1：1–16, 盒子2：17–32, 盒子3：33–48
    %    - 结构：块对角矩阵 [L2L_1, 0, 0; 0, L2L_2, 0; 0, 0, L2L_3]
    %
    % 3. A_layer：
    %    - 行索引：总行数为n_boxes×n_cheb_local，按盒子顺序分块排列
    %      - 盒子k的行：(k-1)×n_cheb_local + 1 到 k×n_cheb_local
    %      - 示例：n_boxes=3, n_cheb_local=16，总行数48
    %        盒子1：1–16, 盒子2：17–32, 盒子3：33–48
    %    - 列索引：总列数为n_boxes×n_neighbors_per_box×n_cheb_local，按盒子顺序分块排列
    %      - 盒子k的列：(k-1)×n_neighbors_per_box×n_cheb_local + 1 到 k×n_neighbors_per_box×n_cheb_local
    %      - 示例：n_boxes=3, n_neighbors_per_box=4, n_cheb_local=16，总列数192
    %        盒子1：1–64, 盒子2：65–128, 盒子3：129–192
    %    - 结构：块对角矩阵 [A_1, 0, 0; 0, A_2, 0; 0, 0, A_3]
    
    if nargin < 2
        use_parallel = 0; % 默认不使用并行计算
    end
    
    % 验证use_parallel参数
    if ~ismember(use_parallel, [0, 1])
        error('use_parallel must be 0 (serial) or 1 (parallel)');
    end
    
    % 如果启用并行计算且并行池未启动，则启动并行池
    if use_parallel == 1 && isempty(gcp('nocreate'))
        parpool;
    end
    
    % 获取总层数
    n_levels = length(mapped_struct_tree);
    
    % 临时存储每层的结果
    temp_struct = mapped_struct_tree;
    
    % 串行遍历每一层
    for lvl = 1:n_levels
        current_nodes = mapped_struct_tree(lvl).nodes;
        n_boxes = length(current_nodes);
        
        % 初始化每层的矩阵列表
        M2M_list = cell(1, n_boxes);
        L2L_list = cell(1, n_boxes);
        A_list = cell(1, n_boxes);
        
        % 根据use_parallel选择并行或串行组装
        if use_parallel == 1
            % 并行收集每层盒子的矩阵
            parfor box_idx = 1:n_boxes
                M2M_list{box_idx} = current_nodes(box_idx).M2M;
                L2L_list{box_idx} = current_nodes(box_idx).L2L;
                A_list{box_idx} = current_nodes(box_idx).A;
            end
        else
            % 串行收集每层盒子的矩阵
            for box_idx = 1:n_boxes
                M2M_list{box_idx} = current_nodes(box_idx).M2M;
                L2L_list{box_idx} = current_nodes(box_idx).L2L;
                A_list{box_idx} = current_nodes(box_idx).A;
            end
        end
        
        % 组装整层矩阵（按块对角形式）
        % 移除空矩阵以避免blkdiag出错
        M2M_valid = M2M_list(~cellfun('isempty', M2M_list));
        L2L_valid = L2L_list(~cellfun('isempty', L2L_list));
        A_valid = A_list(~cellfun('isempty', A_list));
        
        if ~isempty(M2M_valid)
            temp_struct(lvl).M2M_layer = blkdiag(M2M_valid{:});
        else
            temp_struct(lvl).M2M_layer = [];
        end
        
        if ~isempty(L2L_valid)
            temp_struct(lvl).L2L_layer = blkdiag(L2L_valid{:});
        else
            temp_struct(lvl).L2L_layer = [];
        end
        
        if ~isempty(A_valid)
            temp_struct(lvl).A_layer = blkdiag(A_valid{:});
        else
            temp_struct(lvl).A_layer = [];
        end
    end
    
    % 更新mapped_struct_tree
    mapped_struct_tree = temp_struct;
end