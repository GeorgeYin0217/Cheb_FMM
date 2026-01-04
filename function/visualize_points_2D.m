function visualize_points_2D(mapped_struct, target_level)
    % VISUALIZE_LEVEL 可视化指定层级的子盒边界和映射点集
    % 输入:
    %   mapped_struct - 包含映射结果的struct数组
    %   target_level  - 要可视化的层级编号 (例如 1, 2, ...)

    % 检查输入层级是否有效
    n_levels = length(mapped_struct);
    if target_level < 1 || target_level > n_levels
        error('无效的层级编号。请选择 1 到 %d 之间的值。', n_levels);
    end

    % 获取目标层级的节点数据
    level_data = mapped_struct(target_level);
    n_nodes    = length(level_data.nodes);

    % 创建图形窗口
    figure;
    hold on;

    % 初始化坐标范围
    xMinAll = inf; xMaxAll = -inf;
    yMinAll = inf; yMaxAll = -inf;

    % 遍历每个子盒
    for i = 1:n_nodes
        % 提取子盒边界 [x_left, x_right, y_bottom, y_top]
        bnd      = level_data.nodes(i).boundary;
        x_left   = bnd(1);  x_right  = bnd(2);
        y_bottom = bnd(3);  y_top    = bnd(4);

        % 更新全局坐标范围
        xMinAll = min(xMinAll, x_left);
        xMaxAll = max(xMaxAll, x_right);
        yMinAll = min(yMinAll, y_bottom);
        yMaxAll = max(yMaxAll, y_top);

        % 绘制子盒边界 (红色实线)
        rectangle('Position', [x_left, y_bottom, x_right-x_left, y_top-y_bottom], ...
                  'EdgeColor', 'r', 'LineWidth', 2);

        % 提取并绘制映射后的点集
        pts = level_data.nodes(i).points;
        plot(pts(:,1), pts(:,2), 'b.', 'MarkerSize', 3);
    end

    % 设置图形属性并自动调整坐标轴
    axis equal;
    xlim([xMinAll, xMaxAll]);
    ylim([yMinAll, yMaxAll]);
    xlabel('X');
    ylabel('Y');
    title(sprintf('Level %d - Subbox Boundaries and Mapped Points', target_level));
    hold off;
end