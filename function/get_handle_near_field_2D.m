function w_total = get_handle_near_field_2D(q_1, A_matrix_total, h_unaccept_total, indices_total)
    % 构造近场贡献的句柄函数实现
    % 输入：
    %   q_1 - 权重向量
    %   A_matrix_total - 预计算的总核矩阵
    %   h_unaccept_total - 总不可接受点索引向量
    %   indices_total - 总函数点索引向量
    % 输出：
    %   w_total - 近场贡献向量
    
    total_points = length(q_1);
    w_total = zeros(total_points, 1);
    
    % 提取所有权重
    q_selected = q_1(h_unaccept_total);
    
    % 一次性计算所有近场贡献
    w = A_matrix_total * q_selected;
    
    % 将结果放入对应位置
    w_total(indices_total) = w;
end