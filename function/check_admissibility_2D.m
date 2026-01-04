function is_admissible = check_admissibility_2D(node1, node2, eta)
    % 输入
    % node1,node2:子盒1与子盒2的struct结构数组, 含有boundary字段, 表示四边形顶点
    % boundary规则:[左边界x,右边界x,下边界y,上边界y]
    b1 = node1.boundary;
    b2 = node2.boundary;
    
    diam1 = max(b1(2) - b1(1), b1(4) - b1(3));
    diam2 = max(b2(2) - b2(1), b2(4) - b2(3));
    diam = max(diam1, diam2);
    
    dx = max(0, max(b1(1) - b2(2), b2(1) - b1(2)));
    dy = max(0, max(b1(3) - b2(4), b2(3) - b1(4)));
    dist = sqrt(dx^2 + dy^2);
    
    is_admissible = (dist > eta * diam);
    
% % 更详细的调试信息
%     fprintf('node1 boundary: [%f, %f, %f, %f]\n', b1(1), b1(2), b1(3), b1(4));
%     fprintf('node2 boundary: [%f, %f, %f, %f]\n', b2(1), b2(2), b2(3), b2(4));
%     fprintf('diam1: %f, diam2: %f, diam: %f\n', diam1, diam2, diam);
%     fprintf('dx: %f, dy: %f, dist: %f\n', dx, dy, dist);
%     fprintf('eta: %f, eta*diam: %f, is_admissible: %d\n', eta, eta * diam, is_admissible);
end
