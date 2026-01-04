function [tree, level_nodes_struct] = mixquadtree_2D_decomposition(unique_coords, unicoord_indics, max_points, min_points, concentration_threshold,ar_thresh)
    % 输入：
    % unique_coords: 2D点集，[N×2]矩阵，每行是一个点的x,y坐标
    % unicoord_indics: 点集索引集合，[N×1]向量，与unique_coords对应
    % max_points: 每个节点允许的最大点数，超过则继续划分
    % min_points: 每个节点允许的最小点数，避免过度划分
    % concentration_threshold: 集中性判断阈值（例如2/3），为空则不检查
    % 输出：
    % tree: 原始树状结构（含子节点）
    % level_nodes_struct: 重构后的每层节点数组（便于后续处理）
    % ar_thresh:长宽比阈值，设置为很大时不进行二叉树划分

    % ---------------- 参数检查 ----------------
    if size(unique_coords,1) ~= length(unicoord_indics)
        error('unique_coords 和 unicoord_indics 的长度必须一致');
    end
    if ~isempty(concentration_threshold) && (~isnumeric(concentration_threshold) || concentration_threshold <= 0 || concentration_threshold >= 1)
        error('concentration_threshold 必须是(0,1)区间内的数');
    end

    % ---------------- 初始化根节点 ----------------
    x_min = min(unique_coords(:,1)); x_max = max(unique_coords(:,1));
    y_min = min(unique_coords(:,2)); y_max = max(unique_coords(:,2));
    b = [x_min, x_max, y_min, y_max];
    b(1:2) = sort(b(1:2)); b(3:4) = sort(b(3:4));

    level_counts = zeros(1,1); % 每层节点计数器
    level_counts(1) = 1;

    root = struct();
    root.points = unique_coords;
    root.indices = unicoord_indics;
    root.boundary = b;
    root.level_number = [1, 1];
    root.parent_level_number = [0, 0];
    root.leaf_box = 0;

    % ---------------- 开始递归划分 ----------------
    tree = build_node(root);

    % ---------------- 重构成每层数组 ----------------
    level_nodes_struct = reconstruct_tree(tree);

% ---------------- 嵌套函数：递归构建节点 ----------------
function node = build_node(node)
    % 确保边界排序
    bb = node.boundary;
    bb(1:2) = sort(bb(1:2)); 
    bb(3:4) = sort(bb(3:4));
    node.boundary = bb;
    x_min = bb(1); x_max = bb(2);
    y_min = bb(3); y_max = bb(4);
    w = x_max - x_min;
    h = y_max - y_min;
    lvl = node.level_number(1);

    % 叶子节点判断
    if size(node.points,1) <= max_points || size(node.points,1) <= min_points
        node.leaf_box = 1;
        if ~isempty(concentration_threshold)
            node = shrink_boundary(node, concentration_threshold);
        end
        return;
    end

    % 动态中位数划分
    x_mid = median(node.points(:,1));
    y_mid = median(node.points(:,2));

    % 决定使用二叉还是四叉
    % ar_thresh = 2; % 长宽比阈值，可以自己设置
    x = node.points(:,1);
    y = node.points(:,2);

    if w/h > ar_thresh
        % 横向二分 (x方向)
        use_four_split = false;
        idx_left  = (x <= x_mid);  % 左侧（小于等于中位数）
        idx_right = (x >  x_mid);  % 右侧（大于中位数）

        if ~any(idx_left) || ~any(idx_right)
            node.leaf_box = 1;
            if ~isempty(concentration_threshold)
                node = shrink_boundary(node, concentration_threshold);
            end
            return;
        end

        idx_list = {idx_left, idx_right};
        splits = {
            [x_min, x_mid, y_min, y_max];  % 左
            [x_mid, x_max, y_min, y_max]   % 右
        };

    elseif h/w > ar_thresh
        % 纵向二分 (y方向)
        use_four_split = false;
        idx_bottom = (y <= y_mid); % 下侧（小于等于中位数）
        idx_top    = (y >  y_mid); % 上侧（大于中位数）

        if ~any(idx_bottom) || ~any(idx_top)
            node.leaf_box = 1;
            if ~isempty(concentration_threshold)
                node = shrink_boundary(node, concentration_threshold);
            end
            return;
        end

        idx_list = {idx_bottom, idx_top};
        splits = {
            [x_min, x_max, y_min, y_mid];  % 下
            [x_min, x_max, y_mid, y_max]   % 上
        };

    else
        % 四叉划分
        use_four_split = true;
        idx1 = (x <= x_mid) & (y <= y_mid); % 左下
        idx2 = (x >  x_mid) & (y <= y_mid); % 右下
        idx3 = (x >  x_mid) & (y >  y_mid); % 右上
        idx4 = (x <= x_mid) & (y >  y_mid); % 左上

        if ~any(idx1) || ~any(idx2) || ~any(idx3) || ~any(idx4)
            node.leaf_box = 1;
            if ~isempty(concentration_threshold)
                node = shrink_boundary(node, concentration_threshold);
            end
            return;
        end

        idx_list = {idx1, idx2, idx3, idx4};
        splits = {
            [x_min, x_mid, y_min, y_mid];  % 左下
            [x_mid, x_max, y_min, y_mid];  % 右下
            [x_mid, x_max, y_mid, y_max];  % 右上
            [x_min, x_mid, y_mid, y_max]   % 左上
        };
    end

    % 继续生成子节点
    node.children = cell(numel(splits),1);

    % 确保 level_counts 能存储下一层计数
    child_level = lvl + 1;
    if length(level_counts) < child_level
        level_counts(child_level) = 0;
    end

    for i = 1:numel(splits)
        sb = splits{i};
        idx = idx_list{i};
        if ~any(idx)
            continue;
        end

        % 子编号递增
        level_counts(child_level) = level_counts(child_level) + 1;
        child_id = level_counts(child_level);

        % 创建子节点
        child = struct();
        child.points = node.points(idx,:);
        child.indices = node.indices(idx);
        child.boundary = sb;
        child.level_number = [child_level, child_id];
        child.parent_level_number = node.level_number;
        child.leaf_box = 0;

        node.children{i} = build_node(child);
    end
end


    % ---------------- 嵌套函数：边界收缩 ----------------
    function node = shrink_boundary(node, thresh)
        pts = node.points;
        bb = node.boundary;
        bb(1:2) = sort(bb(1:2)); bb(3:4) = sort(bb(3:4));
        x_min = bb(1); x_max = bb(2);
        y_min = bb(3); y_max = bb(4);
        w = x_max - x_min; h = y_max - y_min;
        xp = pts(:,1); yp = pts(:,2);
        cond = all(xp <= x_min + thresh*w) || all(xp >= x_max - thresh*w) || ...
               all(yp <= y_min + thresh*h) || all(yp >= y_max - thresh*h);
        if cond
            nx_min = min(xp); nx_max = max(xp);
            ny_min = min(yp); ny_max = max(yp);
            ext_x = w * 0.05; ext_y = h * 0.05;
            sb = [nx_min-ext_x, nx_max+ext_x, ny_min-ext_y, ny_max+ext_y];
            sb(1:2) = sort(sb(1:2)); sb(3:4) = sort(sb(3:4));
            sb(1) = max(sb(1), x_min);
            sb(2) = min(sb(2), x_max);
            sb(3) = max(sb(3), y_min);
            sb(4) = min(sb(4), y_max);
            node.boundary = sb;
        end
    end

    % ---------------- 嵌套函数：重构层级节点 ----------------
    function level_nodes_struct = reconstruct_tree(tree)
        temp = traverse(tree);
        max_level = max([temp.level]);

        prev_level_map = containers.Map(); % 保存上一层 key->编号 的映射
        level_nodes_struct = repmat(struct('level',[],'nodes',[]), max_level,1);

        for lvl = 1:max_level
            cur = temp([temp.level]==lvl);
            inherited = temp([temp.level]<lvl & [temp.leaf_box]==1);
            nodes = [cur, inherited];

            this_level_map = containers.Map();
            for k = 1:length(nodes)
                nodes(k).new_id = k;
                key = make_key(nodes(k).boundary, nodes(k).indices);
                this_level_map(key) = k;
            end

            N = repmat(struct('level',[],'level_number',[],'parent_level_number',[], ...
                              'boundary',[],'points',[],'indices',[],'leaf_box',[]), 1, length(nodes));
            for k = 1:length(nodes)
                N(k).level = lvl;
                N(k).level_number = nodes(k).new_id;
                if nodes(k).level == lvl
                    N(k).parent_level_number = nodes(k).parent_level_number;
                else
                    prev_key = make_key(nodes(k).boundary, nodes(k).indices);
                    if ~isKey(prev_level_map, prev_key)
                        error('找不到继承节点的父节点，Level %d', lvl);
                    end
                    parent_id = prev_level_map(prev_key);
                    N(k).parent_level_number = [lvl-1, parent_id];
                end
                N(k).boundary = nodes(k).boundary;
                N(k).points = nodes(k).points;
                N(k).indices = nodes(k).indices;
                N(k).leaf_box = nodes(k).leaf_box;
            end
            level_nodes_struct(lvl).level = lvl;
            level_nodes_struct(lvl).nodes = N;

            prev_level_map = this_level_map;
        end

        function key = make_key(boundary, indices)
            boundary = boundary(:)';
            indices = sort(indices(:))';
            key = mat2str([boundary, indices]);
        end

        function L = traverse(node)
            if isempty(node), L = []; return; end
            entry.level = node.level_number(1);
            entry.parent_level_number = node.parent_level_number;
            entry.boundary = node.boundary;
            entry.points = node.points;
            entry.indices = node.indices;
            entry.leaf_box = node.leaf_box;
            L = entry;
            if isfield(node,'children') && ~isempty(node.children)
                for m = 1:length(node.children)
                    L = [L, traverse(node.children{m})];
                end
            end
        end
    end
end