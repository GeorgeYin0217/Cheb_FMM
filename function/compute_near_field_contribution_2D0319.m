%% 后续计算函数
function w_total = compute_near_field_contribution_2D0319(A_matrices_grouped, h_unaccept_grouped, indices_grouped, q_1, use_parallel)
    % 计算FMM近场贡献（自适应分组版本）
    % 输入：
    %   A_matrices_grouped - 分组后的近场矩阵单元数组
    %   h_unaccept_grouped - 分组后的不可接受点索引单元数组
    %   indices_grouped - 分组后的函数点索引单元数组
    %   q_1 - 权重向量
    %   use_parallel - 并行计算开关（0:串行，1:并行）
    % 输出：
    %   w_total - 近场贡献向量
    
    num_groups = length(A_matrices_grouped);
    total_points = length(q_1);
    w_total = zeros(total_points, 1, 'double');
    
    if use_parallel == 1 && ~isempty(gcp('nocreate'))
        % 并行计算
        w_temp = cell(num_groups, 1);
        parfor group_idx = 1:num_groups
            A_group = A_matrices_grouped{group_idx};
            h_unaccept = h_unaccept_grouped{group_idx};
            indices = indices_grouped{group_idx};
            w_temp{group_idx} = struct('w', A_group * q_1(h_unaccept), 'indices', indices);
        end
        
        % 合并结果
        for group_idx = 1:num_groups
            w_data = w_temp{group_idx};
            w_total(w_data.indices) = w_data.w;
        end
    else
        % 串行计算
        for group_idx = 1:num_groups
            A_group = A_matrices_grouped{group_idx};
            h_unaccept = h_unaccept_grouped{group_idx};
            indices = indices_grouped{group_idx};
            w_total(indices) = A_group * q_1(h_unaccept);
        end
    end
end
