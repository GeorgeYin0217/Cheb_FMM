function visualize_admissibility_2D(level_admissibility, target_level, target_number, varargin)
    % 可视化指定层级子盒的可接受性关系
    % 输入：
    %   level_admissibility - analyze_level_admissibility_2D的输出结构
    %   target_level - 目标层级
    %   target_number - 目标子盒编号
    % 可选参数：
    %   'show_numbers' - 是否显示子盒编号 (true/false)
    %   'axis_limits' - 坐标轴范围 [xmin xmax ymin ymax]
    
    % 解析可选参数
    p = inputParser;
    addParameter(p, 'show_numbers', true, @islogical);
    addParameter(p, 'axis_limits', [], @(x)isnumeric(x) && numel(x)==4);
    parse(p, varargin{:});
    
    % 查找目标层级
    level_idx = find([level_admissibility.level] == target_level);
    if isempty(level_idx)
        error('指定层级 %d 不存在', target_level);
    end
    level_data = level_admissibility(level_idx);
    
    % 查找目标子盒
    target_idx = find([level_data.nodes.level_number] == target_number);
    if isempty(target_idx)
        error('在层级 %d 中找不到编号 %d 的子盒', target_level, target_number);
    end
    
    % 获取目标子盒的可接受性关系
    adm_matrix = level_data.admissibility_matrix(target_idx, :);
    
    % 创建新图窗
    figure('Color','white');
    hold on;
    axis equal;
    
    % 设置坐标轴
    if ~isempty(p.Results.axis_limits)
        axis(p.Results.axis_limits);
    else
        % 自动计算坐标范围
        all_bounds = vertcat(level_data.nodes.boundary);
        xlims = [min(all_bounds(:,1)) max(all_bounds(:,2))];
        ylims = [min(all_bounds(:,3)) max(all_bounds(:,4))];
        axis([xlims ylims]);
    end
    xlabel('X'); ylabel('Y');
    title(sprintf('层级 %d - 子盒 %d 的可接受性', target_level, target_number));
    
    % 按层级顺序绘制子盒：先可接受区域 -> 不可接受区域 -> 目标子盒
    for layer = 1:3
        for k = 1:length(level_data.nodes)
            node = level_data.nodes(k);
            bounds = node.boundary;
            
            % 确定绘制条件和颜色
            if layer == 1 && adm_matrix(k) && k ~= target_idx  % 可接受区域（底层）
                color = [0 0.4 0.8];    % 蓝色
                lw = 1.5;
            elseif layer == 2 && ~adm_matrix(k) && k ~= target_idx  % 不可接受区域（中间层）
                color = [0.8 0.2 0.2];  % 红色
                lw = 1.5;
            elseif layer == 3 && k == target_idx  % 目标子盒（顶层）
                color = [0.2 0.8 0.2];  % 绿色
                lw = 2;
            else
                continue;  % 跳过不符合当前层级的子盒
            end
            
            % 绘制边界框
            rectangle('Position', [bounds(1), bounds(3), bounds(2)-bounds(1), bounds(4)-bounds(3)],...
                      'EdgeColor', color, 'LineWidth', lw);
            
            % 显示编号
            if p.Results.show_numbers
                text(mean(bounds(1:2)), mean(bounds(3:4)), num2str(node.level_number),...
                    'HorizontalAlignment', 'center', 'Color', color*0.7);
            end
        end
    end
    
    % 添加图例
    legend_entries = {'目标子盒', '可接受区域', '不可接受区域'};
    legend_colors = [0.2 0.8 0.2; 0 0.4 0.8; 0.8 0.2 0.2];
    
    % 创建图例
    h = zeros(3,1);
    for i = 1:3
        h(i) = plot(NaN, NaN, 's', 'MarkerEdgeColor', legend_colors(i,:),...
            'MarkerFaceColor', legend_colors(i,:), 'MarkerSize', 10);
    end
    legend(h, legend_entries, 'Location', 'bestoutside');
    
    hold off;
end