function L = matrix_laplace(m, n, mu, v)
    %   m - 每组的维数
    %   n - 组数

    % 计算尺度参数 b 从方差
    b = sqrt(v / 2);

    exp1 = exprnd(b, m, n);
    exp2 = exprnd(b, m, n);

    L = mu + exp1 - exp2;
end