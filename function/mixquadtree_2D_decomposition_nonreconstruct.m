function tree = mixquadtree_2D_decomposition_nonreconstruct(unique_coords, unicoord_indics, max_points, min_points, concentration_threshold)
    % 输入：
    % unique_coords: 2D点集，[N×2] 矩阵，每行是一个点的 x,y 坐标
    % unicoord_indics: 点集的索引集合，[N×1] 向量，与 unique_coords 对应
    % max_points: 每个节点最大点数，控制划分深度
    % min_points: 每个网格最小点数，避免空网格
    % concentration_threshold: 集中判断标准（例如 2/3），若为空则不检查叶子节点集中性，可以设置为 [] 表示不检查
    % 输出：
    % tree: 混合四叉/二叉树结构，各节点包含点坐标、索引、层级编号(level, id)、父节点编号(parent_level, parent_id)和叶子标记(leaf_box)

    % 参数校验
    if size(unique_coords,1) ~= length(unicoord_indics)
        error('unique_coords 和 unicoord_indics 的长度必须相同');
    end
    if ~isempty(concentration_threshold) && (~isnumeric(concentration_threshold) || concentration_threshold <= 0 || concentration_threshold >= 1)
        error('concentration_threshold 必须是 0 到 1 之间的数值');
    end

    % 初始边界
    x_min = min(unique_coords(:,1)); x_max = max(unique_coords(:,1));
    y_min = min(unique_coords(:,2)); y_max = max(unique_coords(:,2));
    b = [x_min, x_max, y_min, y_max];

    % 长宽比阈值
    ar_thresh = 1e7;

    % 每层编号计数器
    level_counts = [];

    % 构建根节点
    level_counts(1) = 1;
    root = struct();
    root.points = unique_coords;
    root.indices = unicoord_indics;
    root.boundary = b;
    root.level_number = [1, 1];        % 第1层，第1号
    root.parent_level_number = [0, 0];  % 无父节点
    root.leaf_box = 0;

    % 递归构建混合树
    tree = build_node(root);

    function node = build_node(node)
        % 确保边界排序
        bb = node.boundary;
        bb(1:2)=sort(bb(1:2)); bb(3:4)=sort(bb(3:4));
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
        node.children = {};

        % 四叉或二叉划分决策
        if w/h > ar_thresh
            % 横向二叉划分
            splits = { [x_min, x_mid, y_min, y_max], [x_mid, x_max, y_min, y_max] };
            is_four = false;
        elseif h/w > ar_thresh
            % 纵向二叉划分
            splits = { [x_min, x_max, y_min, y_mid], [x_min, x_max, y_mid, y_max] };
            is_four = false;
        else
            % 标准四叉划分，采用原版 quadtree 的象限逻辑
            is_four = true;
        end

        child_level = lvl + 1;
        if length(level_counts) < child_level
            level_counts(child_level) = 0;
        end

        if ~is_four
            % 二叉或多分支逻辑（保留原逻辑）
            for i = 1:numel(splits)
                sb = splits{i}; sb(1:2)=sort(sb(1:2)); sb(3:4)=sort(sb(3:4));
                idx = node.points(:,1) >= sb(1) & node.points(:,1) <= sb(2) & ...
                      node.points(:,2) >= sb(3) & node.points(:,2) <= sb(4);
                if ~any(idx), continue; end
                level_counts(child_level) = level_counts(child_level) + 1;
                child_id = level_counts(child_level);
                child = struct();
                child.points = node.points(idx,:);
                child.indices = node.indices(idx);
                child.boundary = sb;
                child.level_number = [child_level, child_id];
                child.parent_level_number = node.level_number;
                child.leaf_box = 0;
                node.children{end+1} = build_node(child);
            end
        else
            % 四叉划分，严格按照原代码的 <= 和 > 逻辑
            x = node.points(:,1); y = node.points(:,2);
            idx1 = x <= x_mid & y <= y_mid; % 左下
            idx2 = x >  x_mid & y <= y_mid; % 右下
            idx3 = x >  x_mid & y >  y_mid; % 右上
            idx4 = x <= x_mid & y >  y_mid; % 左上
            % 若任一象限为空，则停止划分
            if ~any(idx1) || ~any(idx2) || ~any(idx3) || ~any(idx4)
                node.leaf_box = 1;
                if ~isempty(concentration_threshold)
                    node = shrink_boundary(node, concentration_threshold);
                end
                return;
            end
            % 定义四个子域边界
            bounds = { [x_min, x_mid, y_min, y_mid], ...
                       [x_mid, x_max, y_min, y_mid], ...
                       [x_mid, x_max, y_mid, y_max], ...
                       [x_min, x_mid, y_mid, y_max] };
            idxs   = { idx1, idx2, idx3, idx4 };
            for i = 1:4
                level_counts(child_level) = level_counts(child_level) + 1;
                child_id = level_counts(child_level);
                sb = bounds{i};
                child = struct();
                child.points = node.points(idxs{i}, :);
                child.indices = node.indices(idxs{i});
                child.boundary = sb;
                child.level_number = [child_level, child_id];
                child.parent_level_number = node.level_number;
                child.leaf_box = 0;
                node.children{i} = build_node(child);
            end
        end
    end

    function node = shrink_boundary(node, thresh)
        pts = node.points;
        bb = node.boundary; bb(1:2)=sort(bb(1:2)); bb(3:4)=sort(bb(3:4));
        x_min = bb(1); x_max = bb(2); y_min = bb(3); y_max = bb(4);
        w = x_max - x_min; h = y_max - y_min;
        xp = pts(:,1); yp = pts(:,2);
        cond = all(xp <= x_min + thresh*w) || all(xp >= x_max - thresh*w) || ...
               all(yp <= y_min + thresh*h) || all(yp >= y_max - thresh*h);
        if cond
            nx_min = min(xp); nx_max = max(xp);
            ny_min = min(yp); ny_max = max(yp);
            ext_x = w * 0.05; ext_y = h * 0.05;
            sb = [nx_min-ext_x, nx_max+ext_x, ny_min-ext_y, ny_max+ext_y];
            sb(1:2)=sort(sb(1:2)); sb(3:4)=sort(sb(3:4));
            sb(1)=max(sb(1), x_min);
            sb(2)=min(sb(2), x_max);
            sb(3)=max(sb(3), y_min);
            sb(4)=min(sb(4), y_max);
            node.boundary = sb;
        end
    end
end
