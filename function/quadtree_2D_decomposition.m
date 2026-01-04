function tree = quadtree_2D_decomposition(unique_coords, unicoord_indics, max_points, min_points, method, concentration_threshold)
    % 输入：
    % unique_coords: 2D点集，[N×2]矩阵，每行是一个点的x,y坐标
    % unicoord_indics: 点集的索引集合，[N×1]向量，与unique_coords对应
    % max_points: 每个节点最大点数，控制划分深度
    % min_points: 每个网格最小点数，避免空网格
    % method: 划分方法，1为方法1（最小点数+空区域检查），2为方法2（动态中位数）
    % concentration_threshold: 集中判断标准（例如2/3），若为空则不检查叶子节点集中性,可以设置为[]表示不检查
    % 输出：
    % tree: 四叉树结构，包含点坐标、索引、层级编号(level, number)、父盒索引(parent_level, parent_number)和叶子标记(leaf_box)

    % 检查输入一致性
    if size(unique_coords, 1) ~= length(unicoord_indics)
        error('unique_coords和unicoord_indics的长度必须相同');
    end
    if ~ismember(method, [1, 2])
        error('method必须为1或2');
    end

    if ~isempty(concentration_threshold) && (~isnumeric(concentration_threshold) || concentration_threshold <= 0 || concentration_threshold >= 1)
        error('concentration_threshold必须是0到1之间的数值');
    end
    
    % 获取初始边界
    x_min = min(unique_coords(:,1));
    x_max = max(unique_coords(:,1));
    y_min = min(unique_coords(:,2));
    y_max = max(unique_coords(:,2));
    
    % 创建初始节点，包含坐标和索引，初始层级编号为(1,1)，无父节点设为[0, 0]，初始leaf_box设为0
    root = struct('points', unique_coords, ...
                  'indices', unicoord_indics, ...
                  'boundary', [x_min, x_max, y_min, y_max], ...
                  'children', [], ...
                  'level_number', [1, 1], ... % 根节点为(1,1)
                  'parent_level_number', [0, 0], ... % 根节点无父节点
                  'leaf_box', 0); % 初始设为非叶子节点
                 
    % 递归构建四叉树并编号
    tree = build_quadtree(root, max_points, min_points, method, concentration_threshold);
    
    % 嵌套函数：递归构建四叉树并赋值层级编号、父盒索引和叶子标记
    function node = build_quadtree(node, max_points, min_points, method, concentration_threshold)
        % 如果点数少于阈值，停止划分，标记为叶子节点
        if size(node.points, 1) <= max_points
            node.leaf_box = 1; % 点数少于max_points，标记为叶子节点
            if ~isempty(concentration_threshold)
                % 如果提供了集中判断标准，检查叶子节点是否点集集中，并更新边界
                node = check_concentration(node, concentration_threshold);
            end
            return;
        end
        
        % 根据方法选择划分策略
        if method == 1
            % 方法1：最小点数阈值 + 空区域检查
            if size(node.points, 1) <= min_points
                node.leaf_box = 1; % 点数少于min_points，无法再分，是叶子节点
                if ~isempty(concentration_threshold)
                    node = check_concentration(node, concentration_threshold);
                end
                return;
            end
            
            % 计算几何中点
            x_mid = (node.boundary(1) + node.boundary(2)) / 2;
            y_mid = (node.boundary(3) + node.boundary(4)) / 2;
        else % method == 2
            % 方法2：动态中位数划分
            x = node.points(:,1);
            y = node.points(:,2);
            x_mid = median(x); % x坐标中位数
            y_mid = median(y); % y坐标中位数
        end
        
        % 将点分为四个子区域
        x = node.points(:,1);
        y = node.points(:,2);
        
        % 四个象限的点集
        idx1 = x <= x_mid & y <= y_mid; % 左下
        idx2 = x > x_mid & y <= y_mid;  % 右下
        idx3 = x > x_mid & y > y_mid;   % 右上
        idx4 = x <= x_mid & y > y_mid;  % 左上
        
        % 检查是否所有子区域都有点，若没有则停止划分，标记为叶子节点
        if ~any(idx1) || ~any(idx2) || ~any(idx3) || ~any(idx4)
            node.leaf_box = 1; % 无法有效划分，是叶子节点
            if ~isempty(concentration_threshold)
                node = check_concentration(node, concentration_threshold);
            end
            return;
        end
        
        % 创建子节点数组
        node.children = cell(4,1);
        
        % 获取当前节点的层数和编号
        current_level = node.level_number(1);
        current_number = node.level_number(2);
        
        % 子区域1：左下
        if any(idx1)
            if current_level == 0
                child_num = 1;
            else
                child_num = (current_number - 1) * 4 + 1;
            end
            child1 = struct('points', node.points(idx1,:), ...
                           'indices', node.indices(idx1), ...
                           'boundary', [node.boundary(1), x_mid, node.boundary(3), y_mid], ...
                           'children', [], ...
                           'level_number', [current_level + 1, child_num], ...
                           'parent_level_number', [current_level, current_number], ...
                           'leaf_box', 0); % 初始设为非叶子节点
            node.children{1} = build_quadtree(child1, max_points, min_points, method, concentration_threshold);
        end
        
        % 子区域2：右下
        if any(idx2)
            if current_level == 0
                child_num = 2;
            else
                child_num = (current_number - 1) * 4 + 2;
            end
            child2 = struct('points', node.points(idx2,:), ...
                           'indices', node.indices(idx2), ...
                           'boundary', [x_mid, node.boundary(2), node.boundary(3), y_mid], ...
                           'children', [], ...
                           'level_number', [current_level + 1, child_num], ...
                           'parent_level_number', [current_level, current_number], ...
                           'leaf_box', 0); % 初始设为非叶子节点
            node.children{2} = build_quadtree(child2, max_points, min_points, method, concentration_threshold);
        end
        
        % 子区域3：右上
        if any(idx3)
            if current_level == 0
                child_num = 3;
            else
                child_num = (current_number - 1) * 4 + 3;
            end
            child3 = struct('points', node.points(idx3,:), ...
                           'indices', node.indices(idx3), ...
                           'boundary', [x_mid, node.boundary(2), y_mid, node.boundary(4)], ...
                           'children', [], ...
                           'level_number', [current_level + 1, child_num], ...
                           'parent_level_number', [current_level, current_number], ...
                           'leaf_box', 0); % 初始设为非叶子节点
            node.children{3} = build_quadtree(child3, max_points, min_points, method, concentration_threshold);
        end
        
        % 子区域4：左上
        if any(idx4)
            if current_level == 0
                child_num = 4;
            else
                child_num = (current_number - 1) * 4 + 4;
            end
            child4 = struct('points', node.points(idx4,:), ...
                           'indices', node.indices(idx4), ...
                           'boundary', [node.boundary(1), x_mid, y_mid, node.boundary(4)], ...
                           'children', [], ...
                           'level_number', [current_level + 1, child_num], ...
                           'parent_level_number', [current_level, current_number], ...
                           'leaf_box', 0); % 初始设为非叶子节点
            node.children{4} = build_quadtree(child4, max_points, min_points, method, concentration_threshold);
        end
        
        % 如果有子节点，则当前节点不是叶子节点，保持leaf_box为0
    end

    % 嵌套函数：检查叶子节点点是否集中在盒子某个区域或一边
    function node = check_concentration(node, thresh)
        if size(node.points, 1) == 0
            return;
        end

        x_min = node.boundary(1);
        x_max = node.boundary(2);
        y_min = node.boundary(3);
        y_max = node.boundary(4);
        box_width = x_max - x_min;
        box_height = y_max - y_min;

        x_points = node.points(:,1);
        y_points = node.points(:,2);

        % 判断是否集中在某一边或某个角落
        left_side   = all(x_points <= x_min + thresh * box_width);
        right_side  = all(x_points >= x_max - thresh * box_width);
        bottom_side = all(y_points <= y_min + thresh * box_height);
        top_side    = all(y_points >= y_max - thresh * box_height);

        % 只要集中在一边或者一个角落
        if left_side || right_side || bottom_side || top_side
            % 重新定义新的边界，尽量紧贴这些点
            x_pmin = min(x_points);
            x_pmax = max(x_points);
            y_pmin = min(y_points);
            y_pmax = max(y_points);

            % 稍微扩展一点边界，避免太紧
            x_extend = box_width * 0.05;
            y_extend = box_height * 0.05;

            new_x_min = max(x_pmin - x_extend, x_min);
            new_x_max = min(x_pmax + x_extend, x_max);
            new_y_min = max(y_pmin - y_extend, y_min);
            new_y_max = min(y_pmax + y_extend, y_max);

            node.boundary = [new_x_min, new_x_max, new_y_min, new_y_max];
        end
    end

end