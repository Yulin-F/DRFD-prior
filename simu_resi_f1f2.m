function [legend_handle,res_list1, res_list2] = simu_resi_f1f2(A, B, C, D, B_f, D_f, B_d, D_d, K, P, f1s, f2s, d, s, x0, x_target)
    n_plot = 400;
    [ny,nu] = size(D);
    [~,nx] = size(A);
    [Gamma,H_us] = matrix_parity(A,B,C,D,s);
    Gamma_orth = matrix_orth(Gamma);
    times = size(f1s, 2)-s;
    
%% simulation f1
    x = [x0,zeros(nx,times+s)];
    u = zeros(nu,times+s);
    y = zeros(ny,times+s);
    f = f1s;
    for k = 1:times+s
        u(:,k) = -K * (x(:,k) - x_target);
        x(:,k+1) = A*x(:,k) + B*u(:,k) +  B_f*f(:,k) + B_d*d(:,k);
        y(:,k) = C*x(:,k) + D*u(:,k) + D_f*f(:,k) + D_d*d(:,k);
    end
    res_list = [];
    for i = s+1 : times+s
        res = Gamma_orth * (vector_augment(y,i,s) -H_us * vector_augment(u,i,s));
        res_norm  = res' * P * res;
        res_list = [res_list , res_norm];
    end
    res_list1 = res_list;

%% simulation f2
    x = [x0,zeros(nx,times+s)];
    u = zeros(nu,times+s);
    y = zeros(ny,times+s);
    f = f2s;
    for k = 1:times+s
        u(:,k) = -K * (x(:,k) - x_target);
        x(:,k+1) = A*x(:,k) + B*u(:,k) +  B_f*f(:,k) + B_d*d(:,k);
        y(:,k) = C*x(:,k) + D*u(:,k) + D_f*f(:,k) + D_d*d(:,k);

    end
    res_list = [];
    for i = s+1 : times+s
        res = Gamma_orth * (vector_augment(y,i,s) -H_us * vector_augment(u,i,s));
        res_norm  = res' * P * res;
        res_list = [res_list , res_norm];
    end
    res_list2 = res_list;
%% plot

    res_plot1 = res_list1( length(res_list1)/2-n_plot+1:length(res_list1)/2+n_plot);
    res_plot2 = res_list2( length(res_list2)/2-n_plot+1:length(res_list2)/2+n_plot);

    plot(res_plot1, 'm-', 'LineWidth', 0.4); % 粉色实线
    hold on 
    plot(res_plot2, 'b-', 'LineWidth', 0.4); % 蓝色实线

    plot([0 n_plot*2],[1 1], 'r--', 'LineWidth', 1);   %红色虚线

    ylim([0 6]);
    xlim([0 800]);

legend_handle=legend({'$J(r)$ for known $f_1$', '$J(r)$ for unknown $f_2$', '$J_{\mathrm{th}}$'}, 'Interpreter', 'latex', 'Location', 'northwest');
ylabel('$J(r)$ and $J_{\mathrm{th}}$', 'Interpreter', 'latex','FontSize',15);
xlabel('Sample number','Interpreter','latex','FontSize',15);
set(gca, 'FontSize', 13);
end

%% vector_augment
function xs = vector_augment(x, k, s)
    xs = [];
    for i = 0:s
        xs = [xs;x(:,k - s + i)];
    end
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
