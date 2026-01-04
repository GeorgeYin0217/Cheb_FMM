function mapped_struct_tree = compute_M2M_L2L_A_each_2D(S, A_func, cheb_nodes, mapped_struct_tree, use_parallel)
    % 输入：
    % S: 句柄函数，接受两个二维坐标向量(a,b)，返回标量（用于M2M）
    % A_func: 句柄函数，接受两个二维坐标向量(x,y)，返回标量（用于系数矩阵A）
    % cheb_nodes: N×2的标准区间[-1,1]×[-1,1]上的切比雪夫节点坐标
    % mapped_struct_tree: 已包含map_points、map_cheb和local_new_admissible_neighbors字段的结构数组
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
        temp_nodes = current_nodes;
        
        % 根据use_parallel选择并行或串行处理每层的盒子
        if use_parallel == 1
            % 并行计算每层的每个盒子
            parfor box_idx = 1:n_boxes
                % 获取当前盒子的map_points和map_cheb
                map_points = current_nodes(box_idx).map_points;
                map_cheb = current_nodes(box_idx).map_cheb;
                neighbors = current_nodes(box_idx).local_new_admissible_neighbors;
                
                % 计算M2M和L2L
                if ~isempty(map_points)
                    n_cheb = size(cheb_nodes, 1);
                    n_map = size(map_points, 1);
                    M2M = zeros(n_cheb, n_map);
                    for i = 1:n_cheb
                        for j = 1:n_map
                            M2M(i,j) = S(cheb_nodes(i,:), map_points(j,:));
                        end
                    end
                    L2L = M2M';
                    temp_nodes(box_idx).M2M = M2M;
                    temp_nodes(box_idx).L2L = L2L;
                else
                    temp_nodes(box_idx).M2M = [];
                    temp_nodes(box_idx).L2L = [];
                end
                
                % 计算系数矩阵A
                if ~isempty(map_cheb) && ~isempty(neighbors)
                    n_cheb_local = size(map_cheb, 1); % 当前盒子的map_cheb点数
                    neighbor_points = [];
                    
                    % 收集所有邻居盒子的map_cheb点
                    for neighbor_idx = neighbors
                        if neighbor_idx >= 1 && neighbor_idx <= n_boxes
                            neighbor_points = [neighbor_points; current_nodes(neighbor_idx).map_cheb];
                        end
                    end
                    
                    if ~isempty(neighbor_points)
                        n_neighbor = size(neighbor_points, 1); % 所有邻居点的总数
                        A_matrix = zeros(n_cheb_local, n_neighbor);
                        
                        % 逐对计算A(x,y)
                        for i = 1:n_cheb_local
                            for j = 1:n_neighbor
                                A_matrix(i,j) = A_func(map_cheb(i,:), neighbor_points(j,:));
                            end
                        end
                        temp_nodes(box_idx).A = A_matrix;
                    else
                        temp_nodes(box_idx).A = [];
                    end
                else
                    temp_nodes(box_idx).A = [];
                end
            end
        else
            % 串行计算每层的每个盒子
            for box_idx = 1:n_boxes
                map_points = current_nodes(box_idx).map_points;
                map_cheb = current_nodes(box_idx).map_cheb;
                neighbors = current_nodes(box_idx).local_new_admissible_neighbors;
                
                % 计算M2M和L2L
                if ~isempty(map_points)
                    n_cheb = size(cheb_nodes, 1);
                    n_map = size(map_points, 1);
                    M2M = zeros(n_cheb, n_map);
                    for i = 1:n_cheb
                        for j = 1:n_map
                            M2M(i,j) = S(cheb_nodes(i,:), map_points(j,:));
                        end
                    end
                    L2L = M2M';
                    temp_nodes(box_idx).M2M = M2M;
                    temp_nodes(box_idx).L2L = L2L;
                else
                    temp_nodes(box_idx).M2M = [];
                    temp_nodes(box_idx).L2L = [];
                end
                
                % 计算系数矩阵A
                if ~isempty(map_cheb) && ~isempty(neighbors)
                    n_cheb_local = size(map_cheb, 1);
                    neighbor_points = [];
                    
                    for neighbor_idx = neighbors
                        if neighbor_idx >= 1 && neighbor_idx <= n_boxes
                            neighbor_points = [neighbor_points; current_nodes(neighbor_idx).map_cheb];
                        end
                    end
                    
                    if ~isempty(neighbor_points)
                        n_neighbor = size(neighbor_points, 1);
                        A_matrix = zeros(n_cheb_local, n_neighbor);
                        
                        for i = 1:n_cheb_local
                            for j = 1:n_neighbor
                                A_matrix(i,j) = A_func(map_cheb(i,:), neighbor_points(j,:));
                            end
                        end
                        temp_nodes(box_idx).A = A_matrix;
                    else
                        temp_nodes(box_idx).A = [];
                    end
                else
                    temp_nodes(box_idx).A = [];
                end
            end
        end
        
        % 更新当前层的nodes
        mapped_struct_tree(lvl).nodes = temp_nodes;
    end
end