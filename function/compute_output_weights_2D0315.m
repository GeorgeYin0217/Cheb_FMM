function fmm_tree = compute_output_weights_2D0315(q_1, mapped_struct_tree, p_end)
    % 输入：
    % q_1: 初始权重向量
    % mapped_struct_tree: 包含每一层信息的结构树
    % p_end: 计算终止层（从上到下计数，1为顶层）
    
    % 初始化输出树结构
    fmm_tree = mapped_struct_tree;
    n_levels = length(mapped_struct_tree);
    
    % 检查输入参数
    if p_end < 1 || p_end > n_levels
        error('p_end must be between 1 and the total number of levels (%d)', n_levels);
    end
    
    % 从最底层（n_levels）开始计算，向上到p_end
    total_m_out_prev = []; % 前一层的total_m_out，用于非叶子计算
    
    for lvl = n_levels:-1:p_end
        current_nodes = mapped_struct_tree(lvl).nodes;
        n_boxes = length(current_nodes);
        
        % 如果是最底层，全是叶子盒
        if lvl == n_levels
            % 提取所有盒子的indices，按顺序拼接q_1对应权重
            q_leaf = zeros(sum(arrayfun(@(node) length(node.indices), current_nodes)), 1);
            offset = 0;
            for box_idx = 1:n_boxes
                indices = current_nodes(box_idx).indices;
                n_indices = length(indices);
                q_leaf(offset + 1 : offset + n_indices) = q_1(indices);
                offset = offset + n_indices;
            end
            
            % 计算total_m_out1 = M2M_leaf * q_leaf
            total_m_out = mapped_struct_tree(lvl).M2M_leaf * q_leaf;
            % 计算total_w_out1 = A_layer * total_m_out1
            total_w_out = mapped_struct_tree(lvl).A_layer * total_m_out;
            
        else
            % 检查当前层的叶子盒和非叶子盒
            leaf_indices = find([current_nodes.leaf_box] == 1); % 叶子盒索引
            n_leaf_boxes = length(leaf_indices);
            
            % 计算m_leaf（如果有叶子盒）
            if ~isempty(mapped_struct_tree(lvl).M2M_leaf) && n_leaf_boxes > 0
                % 提取叶子盒的indices，按顺序拼接q_1对应权重
                q_leaf = zeros(sum(arrayfun(@(idx) length(current_nodes(idx).indices), leaf_indices)), 1);
                offset = 0;
                for i = 1:n_leaf_boxes
                    box_idx = leaf_indices(i);
                    indices = current_nodes(box_idx).indices;
                    n_indices = length(indices);
                    q_leaf(offset + 1 : offset + n_indices) = q_1(indices);
                    offset = offset + n_indices;
                end
                m_leaf = mapped_struct_tree(lvl).M2M_leaf * q_leaf;
            else
                m_leaf = []; % 无叶子盒时为空
            end
            
            % 计算m_nonleaf = M2M_nonleaf * total_m_out_prev
            if ~isempty(total_m_out_prev)
                m_nonleaf = mapped_struct_tree(lvl).M2M_nonleaf * total_m_out_prev;
            else
                m_nonleaf = []; % 前一层无输出时为空
            end
            
            % 拼接total_m_out = [m_leaf; m_nonleaf]
            total_m_out = [m_leaf; m_nonleaf];
            % 计算total_w_out = A_layer * total_m_out
            total_w_out = mapped_struct_tree(lvl).A_layer * total_m_out;
        end
        
        % 存储结果到fmm_tree
        fmm_tree(lvl).total_m_out = total_m_out;
        fmm_tree(lvl).total_w_out = total_w_out;
        
        % 更新total_m_out_prev供下一层使用
        total_m_out_prev = total_m_out;
    end
    
    % 对于未计算的层（1到p_end-1），初始化为空
    for lvl = 1:p_end-1
        fmm_tree(lvl).total_m_out = [];
        fmm_tree(lvl).total_w_out = [];
    end
end