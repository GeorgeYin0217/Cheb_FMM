function level_admissibility = analyze_level_admissibility_2D(tree, eta, use_parallel)
    % 功能：分析四叉树每层的可接受性，并在 nodes 中新增字段存储邻居索引
    % 输入：
    %   tree: 四叉树结构，根节点层级为 1，包含 leaf_box 字段
    %   eta: 可接受性判断参数
    %   use_parallel: 并行计算标识，0 = 不启用，1 = 启用
    % 输出：
    %   level_admissibility: 结构体数组，每层包含：
    %     - level: 层数（从 1 开始）
    %     - nodes: 该层节点数组，新增字段：
    %       - local_new_admissible_neighbors: 局部索引（1,2,3...）
    %       - leaf_box: 叶子盒标识
    %       - h_unaccept: 该盒子的不可接受点索引
    %     - admissibility_matrix: 可接受性矩阵（0=不可接受，1=可接受）

    % 输入参数校验
    if nargin < 3
        use_parallel = 0; % 默认不启用并行计算
    end
    if ~ismember(use_parallel, [0, 1])
        error('use_parallel must be 0 or 1');
    end

    % Step 1: 按层收集所有节点，并获取全局点索引总数
    all_nodes = collect_nodes_by_level(tree);
    max_level = max([all_nodes.level]);
    all_indices = unique(vertcat(all_nodes.indices)); % 全局所有点的索引

    % Step 2: 初始化输出结构体
    level_admissibility = struct('level', {}, 'nodes', {}, 'admissibility_matrix', {});

    % Step 3: 逐层处理，从 1 开始
    for lvl = 1:max_level
        % 获取当前层节点并按 level_number 排序
        level_nodes = all_nodes([all_nodes.level] == lvl);
        if isempty(level_nodes)
            continue;
        end
        [~, sort_idx] = sort([level_nodes.level_number], 'ascend');
        level_nodes = level_nodes(sort_idx);
        num_nodes = length(level_nodes);

        % 获取当前层所有 level_number
        current_level_nums = [level_nodes.level_number];

        % 初始化可接受性矩阵和临时结果存储
        admissibility_matrix = zeros(num_nodes, num_nodes);
        local_new_neighbors_all = cell(num_nodes, 1);
        h_unaccept_all = cell(num_nodes, 1); % 存储不可接受点索引

        if lvl == 1
            % 层级 1：根节点层，无父节点，不可接受点为全部点
            if use_parallel
                parfor i = 1:num_nodes
                    [admiss_row, local_new] = compute_admissibility_lvl1(...
                        i, level_nodes, eta, num_nodes);
                    admissibility_matrix(i, :) = admiss_row;
                    local_new_neighbors_all{i} = local_new;
                    h_unaccept_all{i} = all_indices;
                end
            else
                for i = 1:num_nodes
                    [admiss_row, local_new] = compute_admissibility_lvl1(...
                        i, level_nodes, eta, num_nodes);
                    admissibility_matrix(i, :) = admiss_row;
                    local_new_neighbors_all{i} = local_new;
                    h_unaccept_all{i} = all_indices;
                end
            end
        else
            % 层级 > 1：利用父节点信息
            prev_level_admiss = level_admissibility(lvl - 1);
            prev_admiss_matrix = prev_level_admiss.admissibility_matrix;
            prev_nodes = prev_level_admiss.nodes;
            prev_level_nums = [prev_nodes.level_number];

            % 获取父节点索引
            parent_level_numbers = vertcat(level_nodes.parent_level_number);
            parent_level = parent_level_numbers(:, 1);
            parent_level_num = parent_level_numbers(:, 2);
            parent_indices = zeros(num_nodes, 1);
            for i = 1:num_nodes
                if parent_level(i) == 0
                    parent_indices(i) = 0;
                else
                    idx = find(prev_level_nums == parent_level_num(i), 1);
                    if isempty(idx)
                        error('Level %d, node %d: parent_level_number [%d, %d] not found in level %d', ...
                              lvl, i, parent_level(i), parent_level_num(i), lvl-1);
                    end
                    parent_indices(i) = idx;
                end
            end

            % 计算可接受性矩阵并记录新增邻居和不可接受点索引
            if use_parallel
                parfor i = 1:num_nodes
                    [admiss_row, local_new, h_unaccept] = compute_admissibility_and_h(...
                        i, level_nodes, eta, num_nodes, ...
                        parent_indices, prev_admiss_matrix, prev_nodes);
                    admissibility_matrix(i, :) = admiss_row;
                    local_new_neighbors_all{i} = local_new;
                    h_unaccept_all{i} = h_unaccept;
                end
            else
                for i = 1:num_nodes
                    [admiss_row, local_new, h_unaccept] = compute_admissibility_and_h(...
                        i, level_nodes, eta, num_nodes, ...
                        parent_indices, prev_admiss_matrix, prev_nodes);
                    admissibility_matrix(i, :) = admiss_row;
                    local_new_neighbors_all{i} = local_new;
                    h_unaccept_all{i} = h_unaccept;
                end
            end
        end

        % 将临时结果赋值给 level_nodes
        for i = 1:num_nodes
            level_nodes(i).local_new_admissible_neighbors = local_new_neighbors_all{i};
            level_nodes(i).h_unaccept = h_unaccept_all{i};
        end

        % 保存当前层结果
        level_admissibility(end+1) = struct(...
            'level', lvl, ...
            'nodes', level_nodes, ...
            'admissibility_matrix', admissibility_matrix);
    end
end

% 辅助函数：计算第 1 层的可接受性
function [admiss_row, local_new_neighbors] = compute_admissibility_lvl1(...
    i, level_nodes, eta, num_nodes)
    admiss_row = zeros(1, num_nodes);
    local_new_neighbors = [];
    for j = 1:num_nodes
        if i == j
            admiss_row(j) = 0;
        else
            is_admissible = check_admissibility_2D(level_nodes(i), level_nodes(j), eta);
            admiss_row(j) = is_admissible;
            if is_admissible
                local_new_neighbors(end+1) = j;
            end
        end
    end
end

% 辅助函数：计算第 N 层（N > 1）的可接受性和不可接受点索引
function [admiss_row, local_new_neighbors, h_unaccept] = compute_admissibility_and_h(...
    i, level_nodes, eta, num_nodes, parent_indices, prev_admiss_matrix, prev_nodes)
    admiss_row = zeros(1, num_nodes);
    local_new_neighbors = [];
    parent_i = parent_indices(i);
    
    % 从父节点继承不可接受点索引
    if parent_i == 0
        h_unaccept = level_nodes(i).indices;
    else
        h_unaccept = prev_nodes(parent_i).h_unaccept;
    end

    % 计算可接受性和更新 h_unaccept
    for j = 1:num_nodes
        if i == j
            admiss_row(j) = 0;
        else
            parent_j = parent_indices(j);
            if parent_i == 0 || parent_j == 0
                is_admissible = check_admissibility_2D(level_nodes(i), level_nodes(j), eta);
            elseif prev_admiss_matrix(parent_i, parent_j) == 1
                is_admissible = 1;
            else
                is_admissible = check_admissibility_2D(level_nodes(i), level_nodes(j), eta);
            end
            admiss_row(j) = is_admissible;
            if is_admissible && (parent_i == 0 || prev_admiss_matrix(parent_i, parent_j) == 0)
                local_new_neighbors(end+1) = j;
                h_unaccept = setdiff(h_unaccept, level_nodes(j).indices);
            end
        end
    end
end

% 辅助函数：递归收集所有节点并按层分组，叶子盒递归到最底层并统一编号
function nodes = collect_nodes_by_level(tree)
    nodes = struct('level', {}, 'level_number', {}, 'parent_level_number', {}, ...
                   'boundary', {}, 'points', {}, 'indices', {}, 'leaf_box', {}, 'h_unaccept', {});
    max_level = get_max_level(tree);
    temp_nodes = traverse_tree(tree, []); % 先收集原始节点
    nodes = assign_nodes_to_levels(temp_nodes, max_level); % 再统一分配到每一层

    % 计算树的深度
    function max_lvl = get_max_level(node)
        if isempty(node) || ~isfield(node, 'children') || isempty(node.children)
            max_lvl = node.level_number(1);
            return;
        end
        max_lvl = node.level_number(1);
        for i = 1:length(node.children)
            child_max = get_max_level(node.children{i});
            max_lvl = max(max_lvl, child_max);
        end
    end

    % 递归遍历树，收集原始节点
    function temp_nodes = traverse_tree(node, temp_nodes)
        if isempty(node)
            return;
        end
        if isfield(node, 'leaf_box')
            leaf_box = node.leaf_box;
        else
            leaf_box = [];
        end
        
        temp_nodes(end+1).level = node.level_number(1);
        temp_nodes(end).original_level_number = node.level_number(2); % 保留原始编号
        temp_nodes(end).parent_level_number = node.parent_level_number;
        temp_nodes(end).boundary = node.boundary;
        temp_nodes(end).points = node.points;
        temp_nodes(end).indices = node.indices;
        temp_nodes(end).leaf_box = leaf_box;
        temp_nodes(end).h_unaccept = [];

        if isfield(node, 'children') && ~isempty(node.children)
            for i = 1:length(node.children)
                temp_nodes = traverse_tree(node.children{i}, temp_nodes);
            end
        end
    end

    % 将节点分配到每一层，叶子节点递归到更深层并统一编号
    function nodes = assign_nodes_to_levels(temp_nodes, max_level)
        nodes = struct('level', {}, 'level_number', {}, 'parent_level_number', {}, ...
                       'boundary', {}, 'points', {}, 'indices', {}, 'leaf_box', {}, 'h_unaccept', {});
        level_node_map = cell(max_level, 1); % 映射表：原始标识 -> 新 level_number

        % 第一步：为每一层分配节点和新 level_number
        for lvl = 1:max_level
            % 收集该层的原始节点
            original_nodes = temp_nodes([temp_nodes.level] == lvl);
            % 收集从更浅层继承的叶子节点
            inherited_nodes = struct('level', {}, 'original_level_number', {}, 'parent_level_number', {}, ...
                                     'boundary', {}, 'points', {}, 'indices', {}, 'leaf_box', {}, 'h_unaccept', {});
            for i = 1:length(temp_nodes)
                if temp_nodes(i).leaf_box == 1 && temp_nodes(i).level < lvl
                    inherited_nodes(end+1) = temp_nodes(i);
                end
            end

            % 合并该层节点
            level_nodes = [original_nodes, inherited_nodes];
            num_level_nodes = length(level_nodes);

            % 分配新 level_number 并记录映射
            for i = 1:num_level_nodes
                nodes(end+1) = struct(...
                    'level', lvl, ...
                    'level_number', i, ...
                    'parent_level_number', [], ... % 稍后设置
                    'boundary', level_nodes(i).boundary, ...
                    'points', level_nodes(i).points, ...
                    'indices', level_nodes(i).indices, ...
                    'leaf_box', level_nodes(i).leaf_box, ...
                    'h_unaccept', []);
                % 记录映射
                if lvl == level_nodes(i).level
                    % 原始节点
                    key = [lvl, level_nodes(i).original_level_number];
                else
                    % 继承节点，使用上一层的 key
                    prev_key = [lvl-1, find_prev_level_number(level_nodes(i), lvl-1)];
                    key = prev_key;
                end
                level_node_map{lvl} = [level_node_map{lvl}; {key, i}];
            end
        end

        % 第二步：设置所有节点的 parent_level_number
        for i = 1:length(nodes)
            lvl = nodes(i).level;
            if lvl == 1
                nodes(i).parent_level_number = [0, 0]; % 根节点无父节点
            else
                % 找到该节点在 temp_nodes 中的原始信息
                found = false;
                for j = 1:length(temp_nodes)
                    if isequal(nodes(i).boundary, temp_nodes(j).boundary) && ...
                       isequal(nodes(i).points, temp_nodes(j).points) && ...
                       isequal(nodes(i).indices, temp_nodes(j).indices)
                        if nodes(i).level == temp_nodes(j).level
                            % 原始节点
                            parent_key = temp_nodes(j).parent_level_number;
                            parent_lvl = parent_key(1);
                            parent_num = parent_key(2);
                            % 查找父节点的新编号
                            for k = 1:size(level_node_map{parent_lvl}, 1)
                                if isequal(level_node_map{parent_lvl}{k, 1}, [parent_lvl, parent_num])
                                    nodes(i).parent_level_number = [parent_lvl, level_node_map{parent_lvl}{k, 2}];
                                    found = true;
                                    break;
                                end
                            end
                        else
                            % 继承节点，父节点在上一层
                            prev_lvl = lvl - 1;
                            prev_num = find_prev_level_number(temp_nodes(j), prev_lvl);
                            nodes(i).parent_level_number = [prev_lvl, prev_num];
                            found = true;
                        end
                        break;
                    end
                end
                if ~found
                    error('Cannot set parent_level_number for node at level %d, number %d', ...
                          nodes(i).level, nodes(i).level_number);
                end
            end
        end

        % 辅助函数：找到上一层的 level_number
        function prev_num = find_prev_level_number(node, prev_level)
            for k = 1:length(nodes)
                if nodes(k).level == prev_level && ...
                   isequal(nodes(k).boundary, node.boundary) && ...
                   isequal(nodes(k).points, node.points) && ...
                   isequal(nodes(k).indices, node.indices)
                    prev_num = nodes(k).level_number;
                    return;
                end
            end
            error('Cannot find previous level number for node at level %d', node.level);
        end
    end
end