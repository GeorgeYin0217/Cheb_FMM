function level_admissibility = analyze_level_admissibility_mixtree_2D(level_nodes_struct, eta, use_parallel)
% 功能：分析已重构混合四叉/二叉树每层节点的可接受性，并记录局部新邻居和不可接受点索引。
% 输入：
%   level_nodes_struct: mixquadtree_2D_decomposition输出的每层节点结构体
%   eta: 可接受性控制参数（数值越小要求越严格）
%   use_parallel: 是否启用并行计算（0或1）
% 输出：
%   level_admissibility: 每层节点的可接受性信息，包括：
%      - admissibility_matrix: 节点之间是否可接受的矩阵
%      - nodes: 节点信息（新增字段 local_new_admissible_neighbors 和 h_unaccept）

if nargin < 3
    use_parallel = 0;
end
if ~ismember(use_parallel, [0,1])
    error('use_parallel must be 0 or 1');
end

max_level = numel(level_nodes_struct);
all_nodes = [level_nodes_struct.nodes]; % 展开所有节点
all_indices = unique(vertcat(all_nodes.indices)); % 汇总所有索引

level_admissibility = struct('level', {}, 'nodes', {}, 'admissibility_matrix', {});

for lvl = 1:max_level
    level_nodes = level_nodes_struct(lvl).nodes;
    num_nodes = numel(level_nodes);

    admissibility_matrix = zeros(num_nodes, num_nodes);
    local_new_neighbors_all = cell(num_nodes,1);
    h_unaccept_all = cell(num_nodes,1);

    if lvl == 1
        % 第一层：没有父节点，h_unaccept初始化为全部索引
        if use_parallel
            parfor i = 1:num_nodes
                [row, local_new] = compute_admissibility_lvl1(i, level_nodes, eta, num_nodes);
                admissibility_matrix(i,:) = row;
                local_new_neighbors_all{i} = local_new;
                h_unaccept_all{i} = all_indices;
            end
        else
            for i = 1:num_nodes
                [row, local_new] = compute_admissibility_lvl1(i, level_nodes, eta, num_nodes);
                admissibility_matrix(i,:) = row;
                local_new_neighbors_all{i} = local_new;
                h_unaccept_all{i} = all_indices;
            end
        end
    else
        % 后续层：继承父节点的h_unaccept
        prev_nodes = level_admissibility(lvl-1).nodes;
        prev_admissibility = level_admissibility(lvl-1).admissibility_matrix;

        prev_level_numbers = [prev_nodes.level_number];
        parent_level_numbers = reshape([level_nodes.parent_level_number], 2, []).';
        parent_indices = zeros(num_nodes,1);

        for i = 1:num_nodes
            if all(parent_level_numbers(i,:) == [0,0])
                parent_indices(i) = 0;
            else
                idx = find(prev_level_numbers == parent_level_numbers(i,2), 1);
                if isempty(idx)
                    error('Cannot find parent node at level %d', lvl-1);
                end
                parent_indices(i) = idx;
            end
        end

        if use_parallel
            parfor i = 1:num_nodes
                [row, local_new, h_unaccept] = compute_admissibility_and_h(i, level_nodes, eta, num_nodes, parent_indices, prev_admissibility, prev_nodes);
                admissibility_matrix(i,:) = row;
                local_new_neighbors_all{i} = local_new;
                h_unaccept_all{i} = h_unaccept;
            end
        else
            for i = 1:num_nodes
                [row, local_new, h_unaccept] = compute_admissibility_and_h(i, level_nodes, eta, num_nodes, parent_indices, prev_admissibility, prev_nodes);
                admissibility_matrix(i,:) = row;
                local_new_neighbors_all{i} = local_new;
                h_unaccept_all{i} = h_unaccept;
            end
        end
    end

    % 保存到节点中
    for i = 1:num_nodes
        level_nodes(i).local_new_admissible_neighbors = local_new_neighbors_all{i};
        level_nodes(i).h_unaccept = h_unaccept_all{i};
    end

    % 保存当前层结果
    level_admissibility(lvl) = struct( ...
        'level', lvl, ...
        'nodes', level_nodes, ...
        'admissibility_matrix', admissibility_matrix);
end

end

% ---------------- 子函数：第一层可接受性计算 ----------------
function [admiss_row, local_new_neighbors] = compute_admissibility_lvl1(i, nodes, eta, num_nodes)
    admiss_row = zeros(1, num_nodes);
    local_new_neighbors = [];
    for j = 1:num_nodes
        if i ~= j
            is_adm = check_admissibility_2D(nodes(i), nodes(j), eta);
            admiss_row(j) = is_adm;
            if is_adm
                local_new_neighbors(end+1) = j;
            end
        end
    end
end

% ---------------- 子函数：后续层可接受性计算 ----------------
function [admiss_row, local_new_neighbors, h_unaccept] = compute_admissibility_and_h(i, nodes, eta, num_nodes, parent_indices, prev_matrix, prev_nodes)
    admiss_row = zeros(1, num_nodes);
    local_new_neighbors = [];
    parent_i = parent_indices(i);

    % 继承父节点的不可接受点集合
    if parent_i == 0
        h_unaccept = nodes(i).indices;
    else
        h_unaccept = prev_nodes(parent_i).h_unaccept;
    end

    for j = 1:num_nodes
        if i ~= j
            parent_j = parent_indices(j);
            if parent_i == 0 || parent_j == 0
                is_adm = check_admissibility_2D(nodes(i), nodes(j), eta);
            elseif prev_matrix(parent_i, parent_j) == 1
                is_adm = 1;
            else
                is_adm = check_admissibility_2D(nodes(i), nodes(j), eta);
            end
            admiss_row(j) = is_adm;
            if is_adm && (parent_i == 0 || prev_matrix(parent_i, parent_j) == 0)
                local_new_neighbors(end+1) = j;
                h_unaccept = setdiff(h_unaccept, nodes(j).indices);
            end
        end
    end
end