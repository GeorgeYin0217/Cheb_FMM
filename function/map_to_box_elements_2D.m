function mapped_struct = map_to_box_elements_2D(cheb_nodes, level_admissibility, parallel_flag)
    % MAP_TO_PHYSICAL_ELEMENTS_2D 将标准单元的点映射到层次化矩阵的物理子盒
    % 输入:
    %   cheb_nodes - 标准单元上的点坐标 (n_points × 2)，通常在 [-1, 1] × [-1, 1]
    %   level_admissibility - 层次化矩阵的struct数组
    %   parallel_flag - 并行计算标识 (0:串行, 1:并行)
    % 输出:
    %   mapped_struct - 包含原始数据和新字段 map_cheb 的struct数组

    % 获取层数和点数
    n_levels = length(level_admissibility);
    n_points = size(cheb_nodes, 1);

    % 初始化输出struct数组，复制输入数据
    mapped_struct = level_admissibility;

    % 根据并行标志选择计算方式
    if parallel_flag == 1
        if isempty(gcp('nocreate'))
            parpool; % 启动默认并行池
        end
        % 并行遍历每一层
        parfor lvl = 1:n_levels
            n_nodes = length(level_admissibility(lvl).nodes); % 当前层的子盒数量
            temp_nodes = level_admissibility(lvl).nodes; % 临时存储当前层节点
            for i = 1:n_nodes
                % 提取当前子盒的 boundary
                boundary = temp_nodes(i).boundary; % [x_left, x_right, y_bottom, y_top]
                % 映射到物理子盒
                mapped_points =  map_to_box_2D(cheb_nodes, boundary);
                % 将映射结果添加到新字段 map_cheb
                temp_nodes(i).map_cheb = mapped_points;
            end
            % 更新当前层的 nodes
            mapped_struct(lvl).nodes = temp_nodes;
        end
    else
        % 串行遍历每一层
        for lvl = 1:n_levels
            n_nodes = length(level_admissibility(lvl).nodes); % 当前层的子盒数量
            for i = 1:n_nodes
                % 提取当前子盒的 boundary
                boundary = mapped_struct(lvl).nodes(i).boundary; % [x_left, x_right, y_bottom, y_top]
                % 映射到物理子盒
                mapped_points =  map_to_box_2D(cheb_nodes, boundary);
                % 将映射结果添加到新字段 map_cheb
                mapped_struct(lvl).nodes(i).map_cheb = mapped_points;
            end
        end
    end
end