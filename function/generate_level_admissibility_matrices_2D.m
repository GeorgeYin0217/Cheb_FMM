function admissibility_matrices = generate_level_admissibility_matrices_2D(tree, eta)
% 此为简略版,只输出admissibility_matrices
% 收集各层级节点
    layers = collect_layers(tree);
    
    % 初始化结果矩阵
    max_level = length(layers) - 1;
    admissibility_matrices = cell(max_level + 1, 1);
    
    for level = 0:max_level
        current_layer = layers{level + 1};
        if isempty(current_layer)
            admissibility_matrices{level + 1} = {};
            continue;
        end
        
        % 提取并排序节点编号
        numbers = cellfun(@(n) n.level_number(2), current_layer);
        [~, sorted_idx] = sort(numbers);
        sorted_nodes = current_layer(sorted_idx);
        n = length(sorted_nodes);
        
        % 初始化当前层的邻接矩阵
        level_matrix = cell(n, n);
        
        % 填充矩阵元素
        for i = 1:n
            for j = 1:n
                is_adm = check_admissibility_2D(sorted_nodes{i}, sorted_nodes{j}, eta);
                level_matrix{i,j} = iff(is_adm, '1', '0');
            end
        end
        
        admissibility_matrices{level + 1} = level_matrix;
    end
end

function layers = collect_layers(tree)
    % 广度优先遍历收集各层节点
    layers = {};
    queue = {tree};
    
    while ~isempty(queue)
        node = queue{1};
        queue(1) = [];
        
        % 记录当前层
        level = node.level_number(1);
        if numel(layers) <= level
            layers{level+1} = {};
        end
        layers{level+1}{end+1} = node;
        
        % 子节点入队
        if ~isempty(node.children)
            for c = 1:4
                if ~isempty(node.children{c})
                    queue{end+1} = node.children{c};
                end
            end
        end
    end
end

function str = iff(condition, true_str, false_str)
    % 简化条件判断
    if condition
        str = true_str;
    else
        str = false_str;
    end
end