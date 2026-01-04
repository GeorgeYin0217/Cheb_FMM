function [A_matrices_grouped, h_unaccept_grouped, indices_grouped] = precompute_near_field_matrices_2D0317(tree_6)
    % 预计算FMM近场贡献的分组矩阵和索引向量（自适应分组）
    % 输入：
    %   tree_6 - 包含nodes字段的struct数组（已预计算precomputed_A_matrix）
    % 输出：
    %   A_matrices_grouped - 单元数组，每个元素是分组后的近场矩阵（稀疏格式）
    %   h_unaccept_grouped - 单元数组，每个元素是分组后的不可接受点索引向量
    %   indices_grouped - 单元数组，每个元素是分组后的函数点索引向量
    
    %% 步骤1：获取最底层数据
    num_layers = length(tree_6);
    bottom_layer = tree_6(num_layers).nodes;
    num_boxes = length(bottom_layer); % 注意：应改为length(bottom_layer)
    
    %% 步骤2：估算每个盒子的维度
    rows_per_box = zeros(num_boxes, 1);
    cols_per_box = zeros(num_boxes, 1);
    for box_idx = 1:num_boxes
        box = bottom_layer(box_idx);
        [rows_per_box(box_idx), ~] = size(box.precomputed_A_matrix);
        cols_per_box(box_idx) = length(box.h_unaccept);
    end
    
    % 计算平均维度
    avg_rows = mean(rows_per_box);
    avg_cols = mean(cols_per_box);
    
    %% 步骤3：动态选择分组数量
    % 获取并行核心数
    pool = gcp('nocreate');
    if isempty(pool)
        num_cores = feature('numcores'); % 默认物理核心数
    else
        num_cores = pool.NumWorkers; % 并行池的工作线程数
    end
    
    % 设置内存阈值（例如每个分组矩阵不超过50MB，假设双精度8字节/元素）
    max_bytes_per_group = 50e6; % 50MB
    max_elements_per_group = max_bytes_per_group / 8; % 双精度元素数
    
    % 估算每个分组的目标元素数
    elements_per_box = avg_rows * avg_cols;
    target_boxes_per_group = floor(max_elements_per_group / elements_per_box);
    target_boxes_per_group = max(1, target_boxes_per_group); % 至少1个盒子
    
    % 根据核心数和内存限制确定分组数量
    num_groups = min(num_cores, ceil(num_boxes / target_boxes_per_group));
    group_size = ceil(num_boxes / num_groups); % 每个分组的盒子数
    
    
    %% 步骤4：按组合并矩阵和索引
    A_matrices_grouped = cell(num_groups, 1);
    h_unaccept_grouped = cell(num_groups, 1);
    indices_grouped = cell(num_groups, 1);
    
    for group_idx = 1:num_groups
        % 计算当前组的盒子范围
        start_box = (group_idx - 1) * group_size + 1;
        end_box = min(group_idx * group_size, num_boxes);
        
        % 预估当前组的总行数和列数
        total_rows = sum(rows_per_box(start_box:end_box));
        total_cols = sum(cols_per_box(start_box:end_box));
        
        % 初始化当前组的稀疏矩阵和索引向量
        A_group = sparse(total_rows, total_cols);
        h_unaccept_group = zeros(total_cols, 1);
        indices_group = zeros(total_rows, 1);
        
        % 填充当前组的数据
        row_offset = 0;
        col_offset = 0;
        for box_idx = start_box:end_box
            box = bottom_layer(box_idx);
            A_current = sparse(box.precomputed_A_matrix); % 确保稀疏性
            h_unaccept_current = box.h_unaccept;
            indices_current = box.indices;
            
            [rows, ~] = size(A_current);
            cols = length(h_unaccept_current);
            
            % 填充矩阵和索引
            A_group(row_offset + 1:row_offset + rows, ...
                    col_offset + 1:col_offset + cols) = A_current;
            h_unaccept_group(col_offset + 1:col_offset + cols) = h_unaccept_current;
            indices_group(row_offset + 1:row_offset + rows) = indices_current;
            
            % 更新偏移量
            row_offset = row_offset + rows;
            col_offset = col_offset + cols;
        end
        
        % 存储当前组结果
        A_matrices_grouped{group_idx} = A_group;
        h_unaccept_grouped{group_idx} = h_unaccept_group;
        indices_grouped{group_idx} = indices_group;
    end
end