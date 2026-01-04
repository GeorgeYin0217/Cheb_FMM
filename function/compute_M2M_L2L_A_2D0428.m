function mapped_struct_tree = compute_M2M_L2L_A_2D(S, cheb_nodes, A_func, mapped_struct_tree, bottom_M2M_cols, unique_coords, use_parallel)
    % 输入：
    % S: 句柄函数，接受两个二维坐标向量(a,b)，返回标量，用于M2M
    % cheb_nodes: N×2的标准区间[-1,1]×[-1,1]上的切比雪夫节点坐标
    % A_func: 句柄函数，接受两个二维坐标向量(x,y)，返回标量，用于系数矩阵A
    % mapped_struct_tree: 已包含map_points和map_cheb字段的结构数组
    % bottom_M2M_cols: 底层M2M_layer的列数（map_points总数）
    % unique_coords: 二维坐标点集，用于h_unaccept索引指向的点
    % use_parallel: 0表示串行，1表示并行计算
    
    % 默认参数处理
    if nargin < 7
        use_parallel = 0; % 默认不使用并行计算
    end
    if nargin < 6
        error('unique_coords must be provided as an input parameter');
    end
    if nargin < 5
        error('bottom_M2M_cols must be provided as an input parameter');
    end
      
    % 如果启用并行计算且并行池未启动，则启动并行池
    if use_parallel == 1 && isempty(gcp('nocreate'))
        parpool;
    end
    
    % 获取总层数和切比雪夫节点数
    n_levels = length(mapped_struct_tree);
    n_cheb_per_box = size(cheb_nodes, 1);
    
    % 预定义 child_cols 和 A_col_indices
    child_cols_map = cell(n_levels-1, 1);
    A_col_indices_map = cell(n_levels, 1);
    for lvl = 1:n_levels
        current_nodes = mapped_struct_tree(lvl).nodes;
        A_col_indices_map{lvl} = cell(length(current_nodes), 1);
        for box_idx = 1:length(current_nodes)
            A_col_indices_map{lvl}{box_idx} = compute_A_col_indices(...
                current_nodes(box_idx).local_new_admissible_neighbors, n_cheb_per_box);
        end
        if lvl < n_levels
            next_level_nodes = mapped_struct_tree(lvl+1).nodes;
            child_cols_map{lvl} = cell(length(current_nodes), 1);
            for box_idx = 1:length(current_nodes)
                child_cols = [];
                for child_idx = 1:length(next_level_nodes)
                    if isequal(next_level_nodes(child_idx).parent_level_number, ...
                               [lvl, current_nodes(box_idx).level_number])
                        child_num = next_level_nodes(child_idx).level_number;
                        col_start = (child_num-1) * n_cheb_per_box + 1;
                        col_end = child_num * n_cheb_per_box;
                        child_cols = [child_cols, col_start:col_end];
                    end
                end
                child_cols_map{lvl}{box_idx} = child_cols;
            end
        end
    end
    
    % 主循环：遍历每一层
    for lvl = 1:n_levels
        current_nodes = mapped_struct_tree(lvl).nodes;
        n_boxes = length(current_nodes);
        
        % 计算每层的map_cheb总数
        total_cheb_points = n_boxes * n_cheb_per_box;
        
        % 确定M2M_layer的列数
        if lvl == n_levels
            total_map_points = bottom_M2M_cols;
        else
            total_map_points = length(mapped_struct_tree(lvl+1).nodes) * n_cheb_per_box;
        end
        
        % 预定义矩阵
        M2M_layer = sparse(total_cheb_points, total_map_points);
        A_layer = sparse(total_cheb_points, total_cheb_points);
        temp_nodes = current_nodes;
        
        % 根据use_parallel选择并行或串行处理
        if use_parallel == 1
            temp_results = cell(1, n_boxes);
            parfor box_idx = 1:n_boxes
                temp_results{box_idx} = compute_box_matrices(...
                    current_nodes, box_idx, cheb_nodes, S, A_func, n_cheb_per_box);
                
                % 如果是最底层，计算 precomputed_A_matrix
                if lvl == n_levels
                    box = current_nodes(box_idx);
                    points = box.points;            % 盒子内的函数点坐标 (n_points x 2)
                    h_unaccept = box.h_unaccept;    % 不可接受点索引
                    n_points = size(points, 1);
                    n_unaccept = length(h_unaccept);
                    A_matrix = zeros(n_points, n_unaccept);
                    for i = 1:n_points
                        x = points(i, :);
                        for j = 1:n_unaccept
                            y = unique_coords(h_unaccept(j), :);
                            A_matrix(i, j) = A_func(x, y);
                        end
                    end
                    temp_results{box_idx}.precomputed_A_matrix = A_matrix;
                end
            end
            
            % 处理并行结果
            if lvl == n_levels
                for box_idx = 1:n_boxes
                    temp_struct = temp_results{box_idx};
                    if ~isempty(temp_struct.M2M)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        M2M_layer(row_start:row_end, current_nodes(box_idx).indices) = temp_struct.M2M;
                        temp_nodes(box_idx).M2M = temp_struct.M2M;
                        temp_nodes(box_idx).L2L = temp_struct.L2L;
                    end
                    if ~isempty(temp_struct.A_sub)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        A_layer(row_start:row_end, A_col_indices_map{lvl}{box_idx}) = temp_struct.A_sub;
                        temp_nodes(box_idx).A_sub = temp_struct.A_sub;
                    end
                    % 存储预计算的 A_matrix
                    temp_nodes(box_idx).precomputed_A_matrix = temp_struct.precomputed_A_matrix;
                end
            else
                for box_idx = 1:n_boxes
                    temp_struct = temp_results{box_idx};
                    if ~isempty(temp_struct.M2M)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        M2M_layer(row_start:row_end, child_cols_map{lvl}{box_idx}) = temp_struct.M2M;
                        temp_nodes(box_idx).M2M = temp_struct.M2M;
                        temp_nodes(box_idx).L2L = temp_struct.L2L;
                    end
                    if ~isempty(temp_struct.A_sub)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        A_layer(row_start:row_end, A_col_indices_map{lvl}{box_idx}) = temp_struct.A_sub;
                        temp_nodes(box_idx).A_sub = temp_struct.A_sub;
                    end
                end
            end
            L2L_layer = M2M_layer';
        else
            % 串行处理
            if lvl == n_levels
                for box_idx = 1:n_boxes
                    temp_struct = compute_box_matrices(...
                        current_nodes, box_idx, cheb_nodes, S, A_func, n_cheb_per_box);
                    if ~isempty(temp_struct.M2M)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        M2M_layer(row_start:row_end, current_nodes(box_idx).indices) = temp_struct.M2M;
                        temp_nodes(box_idx).M2M = temp_struct.M2M;
                        temp_nodes(box_idx).L2L = temp_struct.L2L;
                    end
                    if ~isempty(temp_struct.A_sub)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        A_layer(row_start:row_end, A_col_indices_map{lvl}{box_idx}) = temp_struct.A_sub;
                        temp_nodes(box_idx).A_sub = temp_struct.A_sub;
                    end
                    % 计算并存储 precomputed_A_matrix
                    box = current_nodes(box_idx);
                    points = box.points;            % 盒子内的函数点坐标 (n_points x 2)
                    h_unaccept = box.h_unaccept;    % 不可接受点索引
                    n_points = size(points, 1);
                    n_unaccept = length(h_unaccept);
                    A_matrix = zeros(n_points, n_unaccept);
                    for i = 1:n_points
                        x = points(i, :);
                        for j = 1:n_unaccept
                            y = unique_coords(h_unaccept(j), :);
                            A_matrix(i, j) = A_func(x, y);
                        end
                    end
                    temp_nodes(box_idx).precomputed_A_matrix = A_matrix;
                end
            else
                for box_idx = 1:n_boxes
                    temp_struct = compute_box_matrices(...
                        current_nodes, box_idx, cheb_nodes, S, A_func, n_cheb_per_box);
                    if ~isempty(temp_struct.M2M)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        M2M_layer(row_start:row_end, child_cols_map{lvl}{box_idx}) = temp_struct.M2M;
                        temp_nodes(box_idx).M2M = temp_struct.M2M;
                        temp_nodes(box_idx).L2L = temp_struct.L2L;
                    end
                    if ~isempty(temp_struct.A_sub)
                        row_start = (box_idx-1) * n_cheb_per_box + 1;
                        row_end = box_idx * n_cheb_per_box;
                        A_layer(row_start:row_end, A_col_indices_map{lvl}{box_idx}) = temp_struct.A_sub;
                        temp_nodes(box_idx).A_sub = temp_struct.A_sub;
                    end
                end
            end
            L2L_layer = M2M_layer';
        end
        
        % 更新当前层的数据
        mapped_struct_tree(lvl).nodes = temp_nodes;
        mapped_struct_tree(lvl).M2M_layer = M2M_layer;
        mapped_struct_tree(lvl).L2L_layer = L2L_layer;
        mapped_struct_tree(lvl).A_layer = A_layer;
    end
end

% 外部函数：计算单个盒子的M2M、L2L和A_sub（保持不变）
function temp_struct = compute_box_matrices(nodes, box_idx, cheb_nodes, S, A_func, n_cheb_per_box)
    temp_struct = struct();
    
    map_points = nodes(box_idx).map_points;
    if ~isempty(map_points)
        n_map = size(map_points, 1);
        n_cheb = size(cheb_nodes, 1);
        M2M = zeros(n_cheb, n_map);
        for i = 1:n_cheb
            for j = 1:n_map
                M2M(i,j) = S(cheb_nodes(i,:), map_points(j,:));
            end
        end
        L2L = M2M';
        temp_struct.M2M = M2M;
        temp_struct.L2L = L2L;
    else
        temp_struct.M2M = [];
        temp_struct.L2L = [];
    end
    
    x = nodes(box_idx).map_cheb;
    admissible_neighbors = nodes(box_idx).local_new_admissible_neighbors;
    if ~isempty(admissible_neighbors)
        n_neighbors = length(admissible_neighbors);
        y = zeros(n_neighbors * n_cheb_per_box, 2);
        for k = 1:n_neighbors
            neighbor_idx = admissible_neighbors(k);
            start_idx = (k-1) * n_cheb_per_box + 1;
            end_idx = k * n_cheb_per_box;
            y(start_idx:end_idx, :) = nodes(neighbor_idx).map_cheb;
        end
        A_sub = zeros(n_cheb_per_box, n_neighbors * n_cheb_per_box);
        for i = 1:n_cheb_per_box
            for j = 1:(n_neighbors * n_cheb_per_box)
                A_sub(i,j) = A_func(x(i,:), y(j,:));
            end
        end
        temp_struct.A_sub = A_sub;
    else
        temp_struct.A_sub = [];
    end
end

% 辅助函数：计算A_layer的列索引（保持不变）
function col_indices = compute_A_col_indices(admissible_neighbors, n_cheb_per_box)
    if isempty(admissible_neighbors)
        col_indices = [];
        return;
    end
    col_indices = [];
    for k = 1:length(admissible_neighbors)
        neighbor_idx = admissible_neighbors(k);
        col_start = (neighbor_idx-1) * n_cheb_per_box + 1;
        col_end = neighbor_idx * n_cheb_per_box;
        col_indices = [col_indices, col_start:col_end];
    end
end