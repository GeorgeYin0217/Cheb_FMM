function f_final = compute_fmm_2D0317_1(q_1, tree_5, use_parallel, p_end, n_levels)
    % COMPUTE_FMM 执行FMM矩阵-向量乘积计算
    % 输入:
    %   q_1 - 初始向量
    %   tree_5 - 树结构，包含 M2M_layer, A_layer, L2L_layer 等层信息
    %   use_parallel - 并行计算标识 (true/false)
    %   p_end - 预先计算好的中间层索引
    %   n_levels - 树的高度
    % 输出:
    %   f_final - 结果向量，表示一次矩阵乘法的结果

    % 初始化存储中间结果的 cell 数组
    results = cell(n_levels, 3);
    
    % Step 1: 从底层向上计算到 p_end 层
    for lvl = n_levels:-1:p_end
        M_current = tree_5(lvl).M2M_layer;
        A_current = tree_5(lvl).A_layer;
        
        if lvl == n_levels
            m_current = M_current * q_1;
        else
            m_current = M_current * results{lvl+1, 1};
        end
        
        w_current = A_current * m_current;
        
        results{lvl, 1} = m_current; % m 向量
        results{lvl, 2} = w_current; % w 向量
    end
    
    % Step 2: 从 p_end 层向下计算到最后一层
    for lvl = p_end:n_levels
        if lvl == p_end
            l_current = results{lvl, 2};
        else
            L2L_prev = tree_5(lvl-1).L2L_layer;
            l_prev = results{lvl-1, 3};
            w_current = results{lvl, 2};
            l_current = L2L_prev * l_prev + w_current;
        end
        
        results{lvl, 3} = l_current;
    end
    
    % Step 3: 计算远场贡献 f_far
    L2L_last = tree_5(n_levels).L2L_layer;
    l_last = results{n_levels, 3};
    f_far = L2L_last * l_last;
    
    % 计算近场贡献并得到最终结果
    f_near = compute_near_field_contribution_2D0316_1(tree_5, q_1, use_parallel);
    f_final = f_far + f_near;
end