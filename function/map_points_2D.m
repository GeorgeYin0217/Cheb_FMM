function mapped_struct_tree = map_points_2D(mapped_struct_tree, use_parallel)
    % use_parallel: 整数，0表示不启用并行计算，1表示启用并行计算
    if nargin < 2
        use_parallel = 0; % 默认不使用并行计算
    end
     
    % 如果启用并行计算且并行池未启动，则启动并行池
    if use_parallel == 1 && isempty(gcp('nocreate'))
        parpool;
    end
    
    % 获取总层数
    n_levels = length(mapped_struct_tree);
    
    % 遍历每一层
    for lvl = 1:n_levels
        current_nodes = mapped_struct_tree(lvl).nodes;
        n_boxes = length(current_nodes);
        temp_nodes = current_nodes;
        
        % 最底层：直接映射points字段
        if lvl == n_levels
            if use_parallel == 1
                parfor box_idx = 1:n_boxes
                    boundary = current_nodes(box_idx).boundary;
                    points = current_nodes(box_idx).points;
                    
                    x_min = boundary(1);
                    x_max = boundary(2);
                    y_min = boundary(3);
                    y_max = boundary(4);
                    
                    map_points = zeros(size(points));
                    map_points(:,1) = -1 + 2*(points(:,1) - x_min)/(x_max - x_min);
                    map_points(:,2) = -1 + 2*(points(:,2) - y_min)/(y_max - y_min);
                    
                    temp_nodes(box_idx).map_points = map_points;
                end
            else
                for box_idx = 1:n_boxes
                    boundary = current_nodes(box_idx).boundary;
                    points = current_nodes(box_idx).points;
                    
                    x_min = boundary(1);
                    x_max = boundary(2);
                    y_min = boundary(3);
                    y_max = boundary(4);
                    
                    map_points = zeros(size(points));
                    map_points(:,1) = -1 + 2*(points(:,1) - x_min)/(x_max - x_min);
                    map_points(:,2) = -1 + 2*(points(:,2) - y_min)/(y_max - y_min);
                    
                    temp_nodes(box_idx).map_points = map_points;
                end
            end
        else
            % 非最底层：基于子盒子的map_cheb字段映射
            next_level_nodes = mapped_struct_tree(lvl + 1).nodes;
            
            if use_parallel == 1
                parfor box_idx = 1:n_boxes
                    % 获取当前盒子的层级和索引
                    current_level = current_nodes(box_idx).level;
                    current_level_number = current_nodes(box_idx).level_number;
                    boundary = current_nodes(box_idx).boundary;
                    
                    % 提取边界值
                    x_min = boundary(1);
                    x_max = boundary(2);
                    y_min = boundary(3);
                    y_max = boundary(4);
                    
                    % 找到下一层的所有子盒子
                    child_points = [];
                    for child_idx = 1:length(next_level_nodes)
                        parent_info = next_level_nodes(child_idx).parent_level_number;
                        if isequal(parent_info, [current_level, current_level_number])
                            child_points = [child_points; next_level_nodes(child_idx).map_cheb];
                        end
                    end
                    
                    % 映射点集（非底层一定有子盒子，直接使用child_cheb）
                    map_points = zeros(size(child_points));
                    map_points(:,1) = -1 + 2*(child_points(:,1) - x_min)/(x_max - x_min);
                    map_points(:,2) = -1 + 2*(child_points(:,2) - y_min)/(y_max - y_min);
                    
                    temp_nodes(box_idx).map_points = map_points;
                end
            else
                for box_idx = 1:n_boxes
                    current_level = current_nodes(box_idx).level;
                    current_level_number = current_nodes(box_idx).level_number;
                    boundary = current_nodes(box_idx).boundary;
                    
                    x_min = boundary(1);
                    x_max = boundary(2);
                    y_min = boundary(3);
                    y_max = boundary(4);
                    
                    child_points = [];
                    for child_idx = 1:length(next_level_nodes)
                        parent_info = next_level_nodes(child_idx).parent_level_number;
                        if isequal(parent_info, [current_level, current_level_number])
                            child_points = [child_points; next_level_nodes(child_idx).map_cheb];
                        end
                    end
                    
                    % 映射点集（非底层一定有子盒子，直接使用child_cheb）
                    map_points = zeros(size(child_points));
                    map_points(:,1) = -1 + 2*(child_points(:,1) - x_min)/(x_max - x_min);
                    map_points(:,2) = -1 + 2*(child_points(:,2) - y_min)/(y_max - y_min);
                    
                    temp_nodes(box_idx).map_points = map_points;
                end
            end
        end
        
        % 更新当前层的nodes
        mapped_struct_tree(lvl).nodes = temp_nodes;
    end
end