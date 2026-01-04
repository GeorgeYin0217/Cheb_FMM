function quadtree_visualize_2D(tree)
% 绘制2D树结构图像
    figure;
    % 绘制所有点
    plot(tree.points(:,1), tree.points(:,2), 'b.', 'MarkerSize', 10);
    hold on;
    
    % 递归绘制所有节点的边界
    draw_boundaries(tree);
    
    hold off;
    axis equal;
    title('Quadtree Decomposition');
    
    % 嵌套函数：递归绘制边界
    function draw_boundaries(node)
        % 绘制当前节点的边界
        x_min = node.boundary(1);
        x_max = node.boundary(2);
        y_min = node.boundary(3);
        y_max = node.boundary(4);
        plot([x_min, x_max, x_max, x_min, x_min], ...
             [y_min, y_min, y_max, y_max, y_min], 'r-');
        
        % 递归绘制子节点的边界
        if ~isempty(node.children)
            for i = 1:length(node.children)
                if ~isempty(node.children{i})
                    draw_boundaries(node.children{i});
                end
            end
        end
    end
end