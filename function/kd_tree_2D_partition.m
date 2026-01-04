function [clusters, cluster_vertices, all_level_vertices] = kd_tree_2D_partition(points, point_indices, min_size, plot_flag)
    % 输入：
    %   points: 二维点坐标 (N x 2 矩阵)
    %   point_indices: 点索引 (N x 1 向量)
    %   min_size: 最小簇大小 (标量)
    %   plot_flag: 绘图指令 (0 表示不绘图，1 表示绘图)
    % 输出：
    %   clusters: cell 数组，每个 cell 表示一个簇，包含 [x, y, index, cluster_id]
    %   cluster_vertices: cell 数组，每个 cell 表示一个簇的四边形顶点坐标 (4x2 矩阵)
    %   all_level_vertices: cell 数组，存储每一层划分的四边形顶点坐标

    % 检查输入
    if size(points, 2) ~= 2 || length(point_indices) ~= size(points, 1)
        error('输入点集必须是 N x 2 矩阵，且索引长度需匹配点数');
    end

    % 初始化输出变量
    n_points = size(points, 1);
    cluster_ids = zeros(n_points, 1); % 每个点的簇ID
    cluster_vertices = cell(0, 1); % 存储每个簇的顶点
    all_level_vertices = cell(0, 1); % 存储每一层划分的顶点
    current_cluster = 1;

    % 获取初始边界
    parent_bounds = [min(points(:, 1)), max(points(:, 1)), min(points(:, 2)), max(points(:, 2))];

    % 调用递归划分函数
    [cluster_ids, ~, cluster_vertices, all_level_vertices] = partition_space(points, point_indices, min_size, cluster_ids, current_cluster, parent_bounds, cluster_vertices, all_level_vertices);

    % 生成 clusters 输出
    unique_clusters = unique(cluster_ids(cluster_ids > 0));
    clusters = cell(length(unique_clusters), 1);
    for i = 1:length(unique_clusters)
        cluster_idx = cluster_ids == unique_clusters(i);
        clusters{i} = [points(cluster_idx, :), point_indices(cluster_idx), cluster_ids(cluster_idx)];
    end

    % 绘图部分（如果 plot_flag == 1）
    if plot_flag == 1
        figure;
        hold on;

        % 绘制所有点，统一颜色（蓝色），并标注簇序号
        scatter(points(:, 1), points(:, 2), 50, 'b', 'filled'); % 蓝色填充点
        for i = 1:size(points, 1)
            text(points(i, 1), points(i, 2), sprintf('%d', cluster_ids(i)), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                'FontSize', 8, 'Color', 'k');
        end

        % 绘制每一层划分的四边形边界（灰色虚线）
        for i = 1:length(all_level_vertices)
            vertices = all_level_vertices{i};
            plot([vertices(:, 1); vertices(1, 1)], [vertices(:, 2); vertices(1, 2)], 'k--', 'LineWidth', 0.5);
        end

        % 绘制最终簇的四边形边界（黑色实线）
        for i = 1:length(cluster_vertices)
            vertices = cluster_vertices{i};
            plot([vertices(:, 1); vertices(1, 1)], [vertices(:, 2); vertices(1, 2)], 'k-', 'LineWidth', 1);
        end

        title('K-D Tree Space Partitioning with All Level Quadrilaterals');
        xlabel('X');
        ylabel('Y');
        grid on;
        axis equal;
        hold off;
    end

    % 嵌套递归划分函数
    function [cluster_ids, next_cluster, cluster_vertices, all_level_vertices] = partition_space(points, indices, min_size, cluster_ids, current_cluster, parent_bounds, cluster_vertices, all_level_vertices)
        n_points = size(points, 1);
        
        % 记录当前层的四边形顶点
        current_vertices = [parent_bounds(1), parent_bounds(3);  % 左下
                           parent_bounds(2), parent_bounds(3);  % 右下
                           parent_bounds(2), parent_bounds(4);  % 右上
                           parent_bounds(1), parent_bounds(4)]; % 左上
        all_level_vertices{end+1} = current_vertices;

        % 终止条件：点数小于等于 min_size
        if n_points <= min_size
            cluster_ids(indices) = current_cluster;
            cluster_vertices{end+1} = current_vertices;
            next_cluster = current_cluster + 1;
            return;
        end
        
        % 计算当前边界范围
        x_range = max(points(:, 1)) - min(points(:, 1));
        y_range = max(points(:, 2)) - min(points(:, 2));
        
        % 选择分割维度（范围最大的维度）
        if x_range > y_range
            split_dim = 1; % x轴
            [~, sort_idx] = sort(points(:, 1)); % 按x轴排序
        else
            split_dim = 2; % y轴
            [~, sort_idx] = sort(points(:, 2)); % 按y轴排序
        end
        
        % 使用中位数点分割（更符合KD树经典实现）
        mid_idx = floor(n_points / 2) + 1;
        sorted_points = points(sort_idx, :);
        sorted_indices = indices(sort_idx);
        split_value = sorted_points(mid_idx, split_dim);
        
        % 分割点集
        if split_dim == 1
            left_mask = points(:, 1) <= split_value;
            right_mask = points(:, 1) > split_value;
            left_points = points(left_mask, :);
            right_points = points(right_mask, :);
            left_indices = indices(left_mask);
            right_indices = indices(right_mask);
            left_bounds = [parent_bounds(1), split_value, parent_bounds(3), parent_bounds(4)];
            right_bounds = [split_value, parent_bounds(2), parent_bounds(3), parent_bounds(4)];
        else
            left_mask = points(:, 2) <= split_value;
            right_mask = points(:, 2) > split_value;
            left_points = points(left_mask, :);
            right_points = points(right_mask, :);
            left_indices = indices(left_mask);
            right_indices = indices(right_mask);
            left_bounds = [parent_bounds(1), parent_bounds(2), parent_bounds(3), split_value];
            right_bounds = [parent_bounds(1), parent_bounds(2), split_value, parent_bounds(4)];
        end
        
        next_cluster = current_cluster;
        
        % 递归划分左侧
        if ~isempty(left_points)
            [cluster_ids, next_cluster, cluster_vertices, all_level_vertices] = partition_space(left_points, left_indices, min_size, cluster_ids, next_cluster, left_bounds, cluster_vertices, all_level_vertices);
        end
        
        % 递归划分右侧
        if ~isempty(right_points)
            [cluster_ids, next_cluster, cluster_vertices, all_level_vertices] = partition_space(right_points, right_indices, min_size, cluster_ids, next_cluster, right_bounds, cluster_vertices, all_level_vertices);
        end
    end
end