function S = construct_Sab_chebfun_2D(k)
    % 构造四维基函数 S(a, b) 的句柄函数（使用 Chebfun）
    % 输入: k - 切比雪夫插值阶数 (每个方向的点数)
    % 输出: S - 函数句柄，S(a, b) 返回标量值，a和b都是二维向量
    
    % 预计算切比雪夫多项式（Chebfun 对象）
    T = cell(k, 1); % 存储 T_0 到 T_{k-1}
    for n = 0:(k-1)
        T{n+1} = chebpoly(n, [-1, 1]); % 生成第 n 阶切比雪夫多项式
    end
    
    % 定义一维切比雪夫基函数和
    function s = cheb_sum(a, b, k)
        % 输入: a, b - 标量值 (如 a_1 和 b_1)
        % 输出: s - 一维基函数和
        s = T{1}(a) * T{1}(b); % mu = 0 项
        for mu = 1:(k-1)
            s = s + 2 * T{mu+1}(a) * T{mu+1}(b); % mu > 0 项，乘以 2
        end
        s = s / k; % 归一化因子 k^{-1}
    end
    
    % 二维基函数 S(a, b)
    S = @(a, b) cheb_sum(a(1), b(1), k) * cheb_sum(a(2), b(2), k) ;
end