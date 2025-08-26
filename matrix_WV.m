function [W, V] = matrix_WV(A, B, C, D, B_f, D_f, B_d, D_d,s)
    [Gamma,~] = matrix_parity(A,B,C,D,s);
    [~,H_ds] = matrix_parity(A,B_d,C,D_d,s);
    [~,H_fs] = matrix_parity(A,B_f,C,D_f,s);
    Gamma_orth = matrix_orth(Gamma);
    W = Gamma_orth * H_ds;
    V = Gamma_orth * H_fs;

end

%% matrix_parity
function [Gamma_s, Hu_s] = matrix_parity(A, B, C, D, s)
    % 计算扩展可观测性矩阵 H_o,s
    Gamma_s = C;
    for i = 1:s
        Gamma_s = [Gamma_s; C*A^i];
    end
    
    [ny,nu] = size(D);
    
    Hu_s = zeros((s+1)*ny, (s+1)*nu);
    for i = 0:s
        Hu_s(i*ny+1: i*ny+ny, i*nu+1: i*nu+nu) = D;

        for j = 0:s-1-i
            Hu_s((i+j+1)*ny+1: (i+j+2)*ny, i*nu+1: i*nu+nu) = C*A^j*B;
        end
    end
end

%% matrix_orth
function [Gamma_orth] = matrix_orth(Gamma)
    [Q, R] = qr(Gamma);
    r = rank(R);
    Gamma_orth = Q(:, r+1:end);
    Gamma_orth = Gamma_orth';
end
