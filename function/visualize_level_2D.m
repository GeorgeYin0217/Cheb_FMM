function visualize_level_2D(mapped_struct, target_level)
    % VISUALIZE_LEVEL 可视化指定层级的子盒边界和映射点集
    % 输入:
    %   mapped_struct - 包含映射结果的struct数组
    %   target_level - 要可视化的层级编号 (例如  1, 2, ...)
    
    % 检查输入层级是否有效
    n_levels = length(mapped_struct);
    if target_level < 1 || target_level > n_levels
        error('无效的层级编号。请选择 0 到 %d 之间的值。', n_levels );
    end
    
    % 获取目标层级的节点数据
    level_data = mapped_struct(target_level); % MATLAB 索引从 1 开始
    n_nodes = length(level_data.nodes);
    
    % 创建图形窗口
    figure;
    hold on;
    
    % 遍历每个子盒
    for i = 1:n_nodes
        % 提取子盒边界 [x_left, x_right, y_bottom, y_top]
        boundary = level_data.nodes(i).boundary;
        x_left = boundary(1);
        x_right = boundary(2);
        y_bottom = boundary(3);
        y_top = boundary(4);
        
        % 绘制子盒边界 (红色实线)
        rectangle('Position', [x_left, y_bottom, x_right - x_left, y_top - y_bottom], ...
                  'EdgeColor', 'r', 'LineWidth', 2);
        
        % 提取并绘制映射后的点集 (map_cheb)
        mapped_points = level_data.nodes(i).map_cheb;
        plot(mapped_points(:, 1), mapped_points(:, 2), 'b.', 'MarkerSize', 10);
    end
    
    % 设置图形属性
    axis equal; % 保持 x 和 y 轴比例一致
    xlabel('X');
    ylabel('Y');
    title(['Level ', num2str(target_level), ' - Subbox Boundaries and Mapped Points']);
    hold off;
end