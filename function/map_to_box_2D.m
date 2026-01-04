function mapped_points = map_to_box_2D(points, boundary)
    % map_to_box_elements_2D的辅助函数, 将标准单元的点线性映射到物理矩形盒子中
    % 输入:
    %   points - 标准单元上的点坐标 (n_points × 2)，范围 [-1, 1] × [-1, 1]
    %   boundary - 子盒边界 [x_left, x_right, y_bottom, y_top]
    % 输出:
    %   mapped_points - 映射后的物理坐标 (n_points × 2)

    % 提取边界值
    x_left = boundary(1);
    x_right = boundary(2);
    y_bottom = boundary(3);
    y_top = boundary(4);

    % 提取标准坐标
    xi = points(:, 1); % 所有点的 xi 坐标
    eta = points(:, 2); % 所有点的 eta 坐标

    % 向量化线性映射
    mapped_points = zeros(size(points));
    mapped_points(:, 1) = x_left + (x_right - x_left) * (xi + 1) / 2; % x 方向映射
    mapped_points(:, 2) = y_bottom + (y_top - y_bottom) * (eta + 1) / 2; % y 方向映射
end