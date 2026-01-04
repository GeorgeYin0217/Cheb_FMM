function tree = evaluate_fmm_boxquality(tree, coords)
%EVALUATE_FMM_QUALITY 评估FMM树结构的质量指标
%
% 输入参数:
%   tree   - FMM树结构，包含多层节点信息
%   coords - 所有节点的坐标矩阵 (k*2)
%
% 输出参数:
%   tree   - 更新后的树结构，每层增加box_quality字段
%
% box_quality字段包含三个指标向量:
%   - concentration: 点集集中度 (点集包围矩形面积/盒子面积)
%   - aspect_ratio:  盒子长宽比 (长/宽)
%   - max_unaccept_dist: 不可接受点最远距离

    % 遍历树的每一层
    for level = 1:length(tree)
        fprintf('正在处理第%d层，共%d个节点...\n', level, length(tree(level).nodes));
        
        num_nodes = length(tree(level).nodes);
        
        % 初始化三个指标向量
        concentration = zeros(num_nodes, 1);
        aspect_ratio = zeros(num_nodes, 1);
        max_unaccept_dist = zeros(num_nodes, 1);
        
        % 遍历当前层的每个节点
        for node_idx = 1:num_nodes
            node = tree(level).nodes(node_idx);
            
            %% 1. 计算点集集中度
            if ~isempty(node.points)
                % 计算points的包围矩形
                x_coords = node.points(:, 1);
                y_coords = node.points(:, 2);
                
                point_x_min = min(x_coords);
                point_x_max = max(x_coords);
                point_y_min = min(y_coords);
                point_y_max = max(y_coords);
                
                % 点集包围矩形面积
                point_area = (point_x_max - point_x_min) * (point_y_max - point_y_min);
                
                % 处理点集退化的情况（所有点在一条线上或重合）
                if point_area == 0
                    point_area = eps; % 使用极小值避免除零
                end
            else
                point_area = eps; % 空点集情况
            end
            
            % 盒子边界面积
            boundary = node.boundary; % [x_min, x_max, y_min, y_max]
            box_area = (boundary(2) - boundary(1)) * (boundary(4) - boundary(3));
            
            % 点集集中度 = 点集面积 / 盒子面积
            concentration(node_idx) = point_area / box_area;
            
            %% 2. 计算盒子长宽比
            box_width = boundary(2) - boundary(1);
            box_height = boundary(4) - boundary(3);
            
            % 长宽比 = max(宽,高) / min(宽,高)
            aspect_ratio(node_idx) = max(box_width, box_height) / max(min(box_width, box_height), eps);
            
            %% 3. 计算不可接受点最远距离
            if ~isempty(node.h_unaccept) && ~isempty(coords)
                % 获取不可接受点的坐标
                unaccept_coords = coords(node.h_unaccept, :);
                
                % 计算盒子中心
                box_center_x = (boundary(1) + boundary(2)) / 2;
                box_center_y = (boundary(3) + boundary(4)) / 2;
                box_center = [box_center_x, box_center_y];
                
                % 计算所有不可接受点到盒子中心的距离
                distances = sqrt(sum((unaccept_coords - box_center).^2, 2));
                
                % 最远距离
                max_unaccept_dist(node_idx) = max(distances);
            else
                % 没有不可接受点的情况
                max_unaccept_dist(node_idx) = 0;
            end
        end
        
        %% 将计算结果存储到tree结构中
        tree(level).box_quality.concentration = concentration;
        tree(level).box_quality.aspect_ratio = aspect_ratio;
        tree(level).box_quality.max_unaccept_dist = max_unaccept_dist;
        
        % 输出当前层的统计信息
        fprintf('第%d层质量指标统计:\n', level);
        fprintf('  点集集中度: 均值=%.4f, 标准差=%.4f, 范围=[%.4f, %.4f]\n', ...
                mean(concentration), std(concentration), min(concentration), max(concentration));
        fprintf('  盒子长宽比: 均值=%.4f, 标准差=%.4f, 范围=[%.4f, %.4f]\n', ...
                mean(aspect_ratio), std(aspect_ratio), min(aspect_ratio), max(aspect_ratio));
        fprintf('  不可接受点最远距离: 均值=%.4f, 标准差=%.4f, 范围=[%.4f, %.4f]\n\n', ...
                mean(max_unaccept_dist), std(max_unaccept_dist), min(max_unaccept_dist), max(max_unaccept_dist));
    end
    
    fprintf('FMM树质量评估完成！\n');
end