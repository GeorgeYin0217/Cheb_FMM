function f_handle = construct_near_expression_2D(top_level, bottom_level)
    % 构造近场计算的矩阵乘法嵌套表达式
    % 初始化表达式字符串
    expr = '';
    
    % 遍历每个层级，从顶层到底层
    for i = top_level:bottom_level
        % 构造L2L链：从底层到当前层的L2L矩阵乘积
        l2l_chain = '';
        for j = bottom_level:-1:i
            if isempty(l2l_chain)
                l2l_chain = sprintf('tree(%d).L2L_layer', j);
            else
                l2l_chain = sprintf('%s * tree(%d).L2L_layer', l2l_chain, j);
            end
        end
        
        % 构造M2M链：从当前层到底层的M2M矩阵乘积
        m2m_chain = '';
        for j = i:bottom_level
            if isempty(m2m_chain)
                m2m_chain = sprintf('tree(%d).M2M_layer', j);
            else
                m2m_chain = sprintf('%s * tree(%d).M2M_layer', m2m_chain, j);
            end
        end
        m2m_chain = sprintf('%s * q', m2m_chain); % 乘以q
        
        % 当前层贡献项：L2L链 * (A_i * M2M链)
        term = sprintf('(%s) * (tree(%d).A_layer * (%s))', l2l_chain, i, m2m_chain);
        
        % 累加表达式
        if i == top_level
            expr = term;
        else
            expr = sprintf('%s + %s', expr, term);
        end
    end
    
    % 构造匿名函数句柄，输入q并返回矩阵运算结果
    func_str = sprintf('@(q, tree) %s;', expr);
    f_handle = str2func(func_str);
end
