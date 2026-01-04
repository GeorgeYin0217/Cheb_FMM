function S = S_FMM_cheb_2D(a, b, k)
    % 输入：
    % a = [a1, a2]：二维向量 a
    % b = [b1, b2]：二维向量 b
    % k：切比雪夫插值点数（每维）
    % 输出：S(a, b) 的值

    % 提取分量
    a1 = a(1); a2 = a(2);
    b1 = b(1); b2 = b(2);

    % 初始化一维基函数
    S1 = 0; % 第一维 (a1, b1)
    S2 = 0; % 第二维 (a2, b2)

    % 计算一维切比雪夫基函数和
    for mu = 0:(k-1)
        % 使用 Chebfun 的 chebpoly 生成切比雪夫多项式
        T_mu = chebpoly(mu); % 第 mu 阶切比雪夫多项式
        
        % 计算 T_mu 在 a 和 b 上的值
        T_mu_a1 = T_mu(a1);
        T_mu_b1 = T_mu(b1);
        T_mu_a2 = T_mu(a2);
        T_mu_b2 = T_mu(b2);
        
        % 累加：mu=0 系数为 1，mu>0 系数为 2（对称性）
        if mu == 0
            S1 = S1 + T_mu_a1 * T_mu_b1;
            S2 = S2 + T_mu_a2 * T_mu_b2;
        else
            S1 = S1 + 2 * T_mu_a1 * T_mu_b1;
            S2 = S2 + 2 * T_mu_a2 * T_mu_b2;
        end
    end

    % 二维基函数：S(a, b) = S1 * S2 * k^(-2)
    S = S1 * S2 * (k^(-2));
end