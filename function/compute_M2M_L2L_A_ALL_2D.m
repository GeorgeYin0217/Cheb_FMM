function mapped_struct_tree = compute_M2M_L2L_A_ALL_2D(S, cheb_nodes, A_func, mapped_struct_tree, use_parallel)
    % 输入：
    % S: 句柄函数，接受两个二维坐标向量(a,b)，返回标量，用于M2M
    % cheb_nodes: N×2的标准区间[-1,1]×[-1,1]上的切比雪夫节点坐标
    % A_func: 句柄函数，接受两个二维坐标向量(x,y)，返回标量，用于系数矩阵A
    % mapped_struct_tree: 已包含map_points和map_cheb字段的结构数组
    % use_parallel: 0表示串行，1表示并行计算
    
    if nargin < 5
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
    
    % 串行遍历每一层
    for lvl = 1:n_levels
        current_nodes = mapped_struct_tree(lvl).nodes;
        n_boxes = length(current_nodes);
        
        % 计算每层的map_cheb和map_points总数
        n_cheb_per_box = size(current_nodes(1).map_cheb, 1); % 每盒子map_cheb点数，例如25
        total_cheb_points = n_boxes * n_cheb_per_box; % 例如1024
        
        % 计算每盒子map_points点数和总点数
        map_points_counts = zeros(n_boxes, 1);
        for box_idx = 1:n_boxes
            map_points_counts(box_idx) = size(current_nodes(box_idx).map_points, 1);
        end
        total_map_points = sum(map_points_counts); % 所有盒子的map_points总数
        
        % 初始化整层M2M和A矩阵（稀疏矩阵）
        n_cheb = size(cheb_nodes, 1); % 切比雪夫节点数
        M2M_layer = sparse(n_cheb * n_boxes, total_map_points);
        A_layer = sparse(total_cheb_points, total_cheb_points);
        
        % 根据use_parallel选择并行或串行处理
        if use_parallel == 1
            % 创建临时数组存储并行结果
            temp_results = cell(1, n_boxes);
            
            % 并行计算每层的每个盒子
            parfor box_idx = 1:n_boxes
                temp_struct = struct();
                
                % --- 计算单个盒子的M2M ---
                map_points = current_nodes(box_idx).map_points;
                if ~isempty(map_points)
                    n_map = size(map_points, 1);
                    M2M = zeros(n_cheb, n_map);
                    for i = 1:n_cheb
                        for j = 1:n_map
                            M2M(i,j) = S(cheb_nodes(i,:), map_points(j,:));
                        end
                    end
                    temp_struct.M2M = M2M; % 单个盒子的M2M
                else
                    temp_struct.M2M = [];
                end
                
                % --- 计算系数矩阵A的子块 ---
                x = current_nodes(box_idx).map_cheb;
                admissible_neighbors = current_nodes(box_idx).local_new_admissible_neighbors;
                
                if ~isempty(admissible_neighbors)
                    n_neighbors = length(admissible_neighbors);
                    y = zeros(n_neighbors * n_cheb_per_box, 2);
                    y_indices = zeros(1, n_neighbors * n_cheb_per_box);
                    for k = 1:n_neighbors
                        neighbor_idx = admissible_neighbors(k);
                        start_idx = (k-1) * n_cheb_per_box + 1;
                        end_idx = k * n_cheb_per_box;
                        y(start_idx:end_idx, :) = current_nodes(neighbor_idx).map_cheb;
                        neighbor_offset = (neighbor_idx-1) * n_cheb_per_box;
                        y_indices(start_idx:end_idx) = neighbor_offset + 1 : neighbor_offset + n_cheb_per_box;
                    end
                    
                    A1 = zeros(n_cheb_per_box, n_neighbors * n_cheb_per_box);
                    for i = 1:n_cheb_per_box
                        for j = 1:(n_neighbors * n_cheb_per_box)
                            A1(i,j) = A_func(x(i,:), y(j,:));
                        end
                    end
                    
                    x_offset = (box_idx-1) * n_cheb_per_box;
                    x_indices = x_offset + 1 : x_offset + n_cheb_per_box;
                    
                    temp_struct.A1 = A1;
                    temp_struct.x_indices = x_indices;
                    temp_struct.y_indices = y_indices;
                else
                    temp_struct.A1 = [];
                    temp_struct.x_indices = [];
                    temp_struct.y_indices = [];
                end
                
                % 存储到临时数组
                temp_results{box_idx} = temp_struct;
            end
            
            % 并行计算后拼接M2M_layer和A_layer
            map_points_offset = 0;
            for box_idx = 1:n_boxes
                temp_struct = temp_results{box_idx};
                
                % 拼接M2M_layer
                if ~isempty(temp_struct.M2M)
                    n_map = size(temp_struct.M2M, 2);
                    row_start = (box_idx-1) * n_cheb + 1;
                    row_end = box_idx * n_cheb;
                    col_start = map_points_offset + 1;
                    col_end = map_points_offset + n_map;
                    M2M_layer(row_start:row_end, col_start:col_end) = temp_struct.M2M;
                    map_points_offset = map_points_offset + n_map;
                end
                
                % 填充A_layer
                if ~isempty(temp_struct.A1)
                    A_layer(temp_struct.x_indices, temp_struct.y_indices) = temp_struct.A1;
                end
            end
            
            % 直接转置得到L2L_layer
            L2L_layer = M2M_layer';
        else
            % 串行计算每层的每个盒子
            map_points_offset = 0;
            for box_idx = 1:n_boxes
                % --- 计算单个盒子的M2M并拼接 ---
                map_points = current_nodes(box_idx).map_points;
                if ~isempty(map_points)
                    n_map = size(map_points, 1);
                    M2M = zeros(n_cheb, n_map);
                    for i = 1:n_cheb
                        for j = 1:n_map
                            M2M(i,j) = S(cheb_nodes(i,:), map_points(j,:));
                        end
                    end
                    % 直接拼接M2M_layer
                    row_start = (box_idx-1) * n_cheb + 1;
                    row_end = box_idx * n_cheb;
                    col_start = map_points_offset + 1;
                    col_end = map_points_offset + n_map;
                    M2M_layer(row_start:row_end, col_start:col_end) = M2M;
                    map_points_offset = map_points_offset + n_map;
                end
                
                % --- 计算系数矩阵A的子块 ---
                x = current_nodes(box_idx).map_cheb;
                admissible_neighbors = current_nodes(box_idx).local_new_admissible_neighbors;
                
                if ~isempty(admissible_neighbors)
                    n_neighbors = length(admissible_neighbors);
                    y = zeros(n_neighbors * n_cheb_per_box, 2);
                    y_indices = zeros(1, n_neighbors * n_cheb_per_box);
                    for k = 1:n_neighbors
                        neighbor_idx = admissible_neighbors(k);
                        start_idx = (k-1) * n_cheb_per_box + 1;
                        end_idx = k * n_cheb_per_box;
                        y(start_idx:end_idx, :) = current_nodes(neighbor_idx).map_cheb;
                        neighbor_offset = (neighbor_idx-1) * n_cheb_per_box;
                        y_indices(start_idx:end_idx) = neighbor_offset + 1 : neighbor_offset + n_cheb_per_box;
                    end
                    
                    A1 = zeros(n_cheb_per_box, n_neighbors * n_cheb_per_box);
                    for i = 1:n_cheb_per_box
                        for j = 1:(n_neighbors * n_cheb_per_box)
                            A1(i,j) = A_func(x(i,:), y(j,:));
                        end
                    end
                    
                    x_offset = (box_idx-1) * n_cheb_per_box;
                    x_indices = x_offset + 1 : x_offset + n_cheb_per_box;
                    
                    A_layer(x_indices, y_indices) = A1;
                end
            end
            
            % 直接转置得到L2L_layer
            L2L_layer = M2M_layer';
        end
        
        % 更新当前层的数据
        mapped_struct_tree(lvl).nodes = current_nodes; % 不存储M2M和L2L到nodes
        mapped_struct_tree(lvl).M2M_layer = M2M_layer;
        mapped_struct_tree(lvl).L2L_layer = L2L_layer;
        mapped_struct_tree(lvl).A_layer = A_layer;
    end
end