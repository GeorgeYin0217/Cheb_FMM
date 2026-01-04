function mapped_struct_tree = compute_M2M_L2L_A_x_2D(S, cheb_nodes, A_func, mapped_struct_tree, bottom_M2M_cols, unique_coords, use_parallel)
    % 输入：
    % S:    句柄函数，接受四个向量(x1,x2,y1,y2)，返回矩阵，用于M2M/L2L
    % cheb_nodes: N×2 的切比雪夫节点坐标
    % A_func:    句柄函数，接受四个向量(x1,x2,y1,y2)，返回矩阵，用于系数矩阵 A
    % mapped_struct_tree: 已包含 map_points 和 map_cheb 字段的结构数组
    % bottom_M2M_cols:    底层 M2M_layer 的列数
    % unique_coords:      全局坐标点集，用于 precomputed_A_matrix
    % use_parallel:       0=串行，1=并行
    
    if nargin < 7, use_parallel = 0; end
    if nargin < 6, error('unique_coords 必须提供'); end
    if nargin < 5, error('bottom_M2M_cols 必须提供'); end
    
    if use_parallel && isempty(gcp('nocreate'))
        parpool;
    end
    
    n_levels       = length(mapped_struct_tree);
    n_cheb_per_box = size(cheb_nodes,1);
    
    % 预compute child 列映射和 A 列索引映射
    child_cols_map    = cell(n_levels-1,1);
    A_col_indices_map = cell(n_levels,1);
    for lvl = 1:n_levels
        nodes = mapped_struct_tree(lvl).nodes;
        A_col_indices_map{lvl} = cell(length(nodes),1);
        for b = 1:length(nodes)
            A_col_indices_map{lvl}{b} = compute_A_col_indices( ...
                nodes(b).local_new_admissible_neighbors, n_cheb_per_box);
        end
        if lvl < n_levels
            next_nodes = mapped_struct_tree(lvl+1).nodes;
            child_cols_map{lvl} = cell(length(nodes),1);
            for b = 1:length(nodes)
                cols = [];
                for c = 1:length(next_nodes)
                    if isequal(next_nodes(c).parent_level_number, [lvl,nodes(b).level_number])
                        num = next_nodes(c).level_number;
                        cols = [cols, (num-1)*n_cheb_per_box+1 : num*n_cheb_per_box];
                    end
                end
                child_cols_map{lvl}{b} = cols;
            end
        end
    end
    
    % 主循环：每层
    for lvl = 1:n_levels
        nodes     = mapped_struct_tree(lvl).nodes;
        n_boxes   = length(nodes);
        total_cheb = n_boxes * n_cheb_per_box;
        if lvl == n_levels
            total_map = bottom_M2M_cols;
        else
            total_map = length(mapped_struct_tree(lvl+1).nodes) * n_cheb_per_box;
        end
        
        % —— 第一遍：预估非零元总数
        nnzM = 0;
        nnzA = 0;
        for b = 1:n_boxes
            % M2M: map_points 非空时计算
            n_map = size(nodes(b).map_points,1);
            nnzM = nnzM + n_cheb_per_box * n_map;
            % A_sub: 有邻居时
            n_nb = length(nodes(b).local_new_admissible_neighbors);
            nnzA = nnzA + n_cheb_per_box * n_cheb_per_box * n_nb;
        end
        
        % —— 预分配三元组数组
        mi = zeros(nnzM,1); mj = zeros(nnzM,1); mv = zeros(nnzM,1);
        ai = zeros(nnzA,1); aj = zeros(nnzA,1); av = zeros(nnzA,1);
        ptrM = 1;  ptrA = 1;
        
        temp_nodes = nodes;
        
        if use_parallel
            % 并行：先计算各盒矩阵
            temp_results = cell(1,n_boxes);
            parfor b = 1:n_boxes
                temp_results{b} = compute_box_matrices( ...
                    nodes, b, cheb_nodes, S, A_func, n_cheb_per_box );
                if lvl == n_levels
                    box = nodes(b);
                    pts = box.points;
                    h_un = box.h_unaccept;
                    if ~isempty(h_un)
                        px1 = pts(:,1); px2 = pts(:,2);
                        uy1 = unique_coords(h_un,1)'; uy2 = unique_coords(h_un,2)';
                        temp_results{b}.precomputed_A_matrix = A_func(px1,px2,uy1,uy2);
                    else
                        temp_results{b}.precomputed_A_matrix = zeros(size(pts,1),0);
                    end
                end
            end
            
            % 串行聚合三元组
            for b = 1:n_boxes
                tr = temp_results{b};
                rs = (b-1)*n_cheb_per_box + (1:n_cheb_per_box);
                
                % M2M/L2L 三元组
                if ~isempty(tr.M2M)
                    if lvl == n_levels
                        colsM = nodes(b).indices;
                    else
                        colsM = child_cols_map{lvl}{b};
                    end
                    [R,C] = ndgrid(rs, colsM);
                    n = numel(R);
                    mi(ptrM:ptrM+n-1) = R(:);
                    mj(ptrM:ptrM+n-1) = C(:);
                    mv(ptrM:ptrM+n-1) = tr.M2M(:);
                    ptrM = ptrM + n;
                    temp_nodes(b).M2M = tr.M2M;
                    temp_nodes(b).L2L = tr.L2L;
                end
                
                % A_sub 三元组
                if ~isempty(tr.A_sub)
                    colsA = A_col_indices_map{lvl}{b};
                    [R,C] = ndgrid(rs, colsA);
                    n = numel(R);
                    ai(ptrA:ptrA+n-1) = R(:);
                    aj(ptrA:ptrA+n-1) = C(:);
                    av(ptrA:ptrA+n-1) = tr.A_sub(:);
                    ptrA = ptrA + n;
                    temp_nodes(b).A_sub = tr.A_sub;
                end
                
                % precomputed_A_matrix
                if lvl == n_levels
                    temp_nodes(b).precomputed_A_matrix = tr.precomputed_A_matrix;
                end
            end
            
        else
            % 串行：逐盒计算并聚合
            for b = 1:n_boxes
                tr = compute_box_matrices( ...
                    nodes, b, cheb_nodes, S, A_func, n_cheb_per_box );
                rs = (b-1)*n_cheb_per_box + (1:n_cheb_per_box);
                
                % M2M/L2L
                if ~isempty(tr.M2M)
                    if lvl == n_levels
                        colsM = nodes(b).indices;
                    else
                        colsM = child_cols_map{lvl}{b};
                    end
                    [R,C] = ndgrid(rs, colsM);
                    n = numel(R);
                    mi(ptrM:ptrM+n-1) = R(:);
                    mj(ptrM:ptrM+n-1) = C(:);
                    mv(ptrM:ptrM+n-1) = tr.M2M(:);
                    ptrM = ptrM + n;
                    temp_nodes(b).M2M = tr.M2M;
                    temp_nodes(b).L2L = tr.L2L;
                end
                
                % A_sub
                if ~isempty(tr.A_sub)
                    colsA = A_col_indices_map{lvl}{b};
                    [R,C] = ndgrid(rs, colsA);
                    n = numel(R);
                    ai(ptrA:ptrA+n-1) = R(:);
                    aj(ptrA:ptrA+n-1) = C(:);
                    av(ptrA:ptrA+n-1) = tr.A_sub(:);
                    ptrA = ptrA + n;
                    temp_nodes(b).A_sub = tr.A_sub;
                end
                
                % precomputed_A_matrix
                if lvl == n_levels
                    box = nodes(b);
                    pts = box.points;
                    h_un = box.h_unaccept;
                    if ~isempty(h_un)
                        px1 = pts(:,1); px2 = pts(:,2);
                        uy1 = unique_coords(h_un,1)'; uy2 = unique_coords(h_un,2)';
                        temp_nodes(b).precomputed_A_matrix = A_func(px1,px2,uy1,uy2);
                    else
                        temp_nodes(b).precomputed_A_matrix = zeros(size(pts,1),0);
                    end
                end
            end
        end
        
        % 一次性构造稀疏矩阵
        M2M_layer = sparse(mi, mj, mv, total_cheb, total_map);
        A_layer   = sparse(ai, aj, av, total_cheb, total_cheb);
        L2L_layer = M2M_layer';
        
        % 更新结构
        mapped_struct_tree(lvl).nodes     = temp_nodes;
        mapped_struct_tree(lvl).M2M_layer = M2M_layer;
        mapped_struct_tree(lvl).L2L_layer = L2L_layer;
        mapped_struct_tree(lvl).A_layer   = A_layer;
    end
end

% -------------------------------------------------------------------------
function temp = compute_box_matrices(nodes, idx, cheb_nodes, S, A_func, n_cheb)
    % 向量化计算单盒的 M2M/L2L/A_sub
    temp = struct('M2M',[],'L2L',[],'A_sub',[]);
    
    % M2M/L2L
    mp = nodes(idx).map_points;
    if ~isempty(mp)
        x1 = cheb_nodes(:,1); x2 = cheb_nodes(:,2);
        y1 = mp(:,1)';      y2 = mp(:,2)';
        M2M = S(x1,x2,y1,y2);
        temp.M2M = M2M;
        temp.L2L = M2M';
    end
    
    % A_sub
    nbrs = nodes(idx).local_new_admissible_neighbors;
    if ~isempty(nbrs)
        x  = nodes(idx).map_cheb;
        x1 = x(:,1); x2 = x(:,2);
        Y  = zeros(length(nbrs)*n_cheb,2);
        for k = 1:length(nbrs)
            Y((k-1)*n_cheb+1:k*n_cheb,:) = nodes(nbrs(k)).map_cheb;
        end
        y1 = Y(:,1)'; y2 = Y(:,2)';
        temp.A_sub = A_func(x1,x2,y1,y2);
    end
end

% -------------------------------------------------------------------------
function cols = compute_A_col_indices(adm, n_cheb)
    if isempty(adm)
        cols = [];
    else
        cols = arrayfun(@(k) (adm(k)-1)*n_cheb + (1:n_cheb), 1:length(adm), 'uni', false);
        cols = [cols{:}];
    end
end
