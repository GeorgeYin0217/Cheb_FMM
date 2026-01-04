function result = process_element_parallel_2D(i, elements, coordinates, points, weights, n_points)
    % PROCESS_ELEMENT 处理单个单元的映射计算（并行模式）
    % 输入:
    %   i - 当前单元索引
    %   elements - 单元连接信息 (4 × n_elements)
    %   coordinates - 所有节点的物理坐标 (n_nodes × 2)
    %   points - 标准单元上的点坐标 (n_points × 2)
    %   weights - 对应的权重 (n_points × 1)
    %   n_points - 点的数量
    % 输出:
    %   result - 当前单元的投影坐标和相关信息 (n_points × 4)

    element_nodes = elements(:, i);
    element_coordinates = coordinates(element_nodes, :);
    x = element_coordinates(:, 1);
    y = element_coordinates(:, 2);
    mapped_points = zeros(n_points, 2);
    
    for j = 1:n_points
        xi = points(j, 1);
        eta = points(j, 2);
        N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
        mapped_points(j, :) = [N*x, N*y];
    end
    
    J_p = get_2D_jacobian(points, mapped_points);
    result = [mapped_points(:, 1), mapped_points(:, 2), weights, J_p];
end