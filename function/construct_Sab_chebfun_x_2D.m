function S = construct_Sab_chebfun_x_2D(k)
    % 构造向量化的二维切比雪夫基函数 S(a,b)
    % 输入: k - 切比雪夫插值阶数
    % 输出: S - 函数句柄，接受四个矩阵 A1,A2,B1,B2，返回同尺寸矩阵

    % 预计算切比雪夫多项式（Chebfun 对象）
    T = cell(k,1);
    for n = 0:(k-1)
        T{n+1} = chebpoly(n, [-1,1]);
    end

    % 定义向量化的一维基函数和
    function S1 = cheb_sum_vec(A, B)
        % A, B: 同尺寸矩阵，元素为标量坐标
        % 输出 S1: 同尺寸矩阵，元素为基函数和
        % 初始 mu=0 项
        S1 = T{1}(A) .* T{1}(B);
        for mu = 1:(k-1)
            S1 = S1 + 2 .* ( T{mu+1}(A) .* T{mu+1}(B) );
        end
        S1 = S1 ./ k;
    end

    % 二维基函数 S(a1,a2,b1,b2)
    S = @(A1,A2,B1,B2) cheb_sum_vec(A1, B1) .* cheb_sum_vec(A2, B2);
end
