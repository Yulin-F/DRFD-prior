clear;clc;close all;
rng(818);
omega = [0.01 0.1 1 ];
s = 4;  %parity space order
MonteCarlo = 250;
N_single = 800;
times = MonteCarlo*N_single;
alpha = 0.01;
beta = [0.05 0.03 0.01];
sample = 200;
theta_list = 0.001:0.001:0.01;    % Wasserstein radius grid
fold = 5;
var = 0.4;

A=[0.8945 0.0048 0.1005;
   0.0048 0.8500 0.0801;
   0.1005 0.0801 0.8164];
B=[0.0317 0;
   0 0.0309;
   0.0018 0.001];
C=[1 0 0;
   0 0 1];
D=[0 0;
   0 0];
B_d= eye(3);
D_d= [0 0 0;
      0 0 0];
B_f=[1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0];
D_f=[0 0 0 1 0;
     0 0 0 0 1];
x0 = [15 10 13]';
x_target = [30 17 15]';
f1 = [-2 0 0 0 0]';
f2 = [0 5 10 0 0]';

dist = matrix_laplace(size(B_d,2),times+s, 0, var);
[W, V] = matrix_WV(A, B, C, D, B_f, D_f, B_d, D_d, s);
xi_hat = W * matrix_laplace(size(W,2), sample, 0, var);

theta1=0.008;
theta2=0.009;

%% Wasserstein radius calibration 
FAR_total= {};
for rho =1:2
FAR_list = [];  
for theta = theta_list
    n_fold = sample/fold;
    res_list = [];
    for i = 1:fold
        xi_train = xi_hat;
        xi_test = xi_train(:,(i-1)*n_fold+1:i*n_fold);
        xi_train(:,(i-1)*n_fold+1:i*n_fold) = [];
        P = FD_DRO(A, B, C, D, B_f, D_f, B_d, D_d, s, xi_train, theta, alpha ,rho);
        for j = 1:n_fold
            res_list = [res_list; xi_test(:,j)'*P*xi_test(:,j)];
        end
    end
    FAR_list =[FAR_list sum(res_list > 1)/sample];
end
FAR_total{rho} = FAR_list;
end


figure(1)
plot(theta_list,FAR_total{1},'*-','LineWidth', 1);
hold on;
plot(theta_list,FAR_total{2},'o-','LineWidth', 1);

xlabel('$\theta$', 'Interpreter', 'latex','FontSize', 15);
ylabel('Cross-validation FAR', 'Interpreter', 'latex','FontSize', 15);
legend({'DR$1$','DR$2$'}, 'Interpreter', 'latex','FontSize', 15);
grid on
set(1, 'Position', [100, 100, 550, 250]); 

theta1 = min(theta_list(FAR_total{1}<=alpha));
theta2 = min(theta_list(FAR_total{2}<=alpha));


%% FD design
[P10] = FD_initial(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta1, alpha);
%DR1
[P11] = FD_DRO(A, B, C, D, B_f, D_f, B_d, D_d, s, xi_hat, theta1, alpha, 1);
[P12 ,kappa12,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta1, alpha, beta(1), omega(1), 1);
[P13 ,kappa13,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta1, alpha, beta(1), omega(2), 1);
[P14 ,kappa14,obj_list1] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta1, alpha, beta(1), omega(3), 1);
[P15 ,kappa15,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta1, alpha, beta(2), omega(3), 1);
[P16 ,kappa16,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta1, alpha, beta(3), omega(3), 1);

%DR2
[P20] = FD_initial(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta2, alpha);
[P21] = FD_DRO(A, B, C, D, B_f, D_f, B_d, D_d, s, xi_hat, theta2, alpha, 2);
[P22 ,kappa22,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta2, alpha, beta(1), omega(1), 2);
[P23 ,kappa23,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta2, alpha, beta(1), omega(2), 2);
[P24 ,kappa24,obj_list2] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta2, alpha, beta(1), omega(3), 2);
[P25 ,kappa25,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta2, alpha, beta(2), omega(3), 2);
[P26 ,kappa26,~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta2, alpha, beta(3), omega(3), 2);

% Convergence fig
figure(2)
yyaxis left
plot(-obj_list1,'*-','LineWidth', 1);
hold on;
ylabel('$\rho_1(P)+\gamma \eta$', 'Interpreter', 'latex','FontSize', 15);
yyaxis right
plot(-obj_list2,'o-','LineWidth', 1);

ylabel('$\rho_2(P)+\gamma \eta$', 'Interpreter', 'latex','FontSize', 15);
xlabel('Iteration', 'Interpreter', 'latex','FontSize', 15);
legend({'NDR$1$-C','NDR$2$-C'}, 'Interpreter', 'latex','FontSize', 15,'Location','southeast');
grid on
set(gca, 'XTick', 1:1:6)
set(2, 'Position', [100, 100, 560, 210]);  


% Residual fig
figure(3);
f1s = zeros(size(f1,1), times + s); 
f1s(:, times/2 + 1 + s:end) = repmat(f1, 1, times/2); 
f2s = zeros(size(f1,1), times + s); 
f2s(:, times/2 + 1 + s:end) = repmat(f2, 1, times/2); 
K = place(A, B, [0.5 0.3 -0.3]);

subplot_titles = {'Init$1$', 'DR$1$','NDR$1$-A','NDR$1$-B','NDR$1$-C','NDR$1$-D', 'NDR$1$-E', 'Init$2$', 'DR$2$',  'NDR$2$-A',  'NDR$2$-B',  'NDR$2$-C', 'NDR$2$-D', 'NDR$2$-E'};
P = {P10, P11, P12, P13, P14, P15, P16, P20, P21, P22, P23, P24, P25, P26};

FAR1 = zeros(2, 7);
FDR1 = zeros(2, 7);
FAR2 = zeros(2, 7);
FDR2 = zeros(2, 7);

ha = tight_subplot(2, 4, [0.2 0.05], [0.12 0.07], [0.05 0.05]);
set(3, 'Position', [100, 100, 1200, 450]);
for i = 1:7
    axes(ha(i)); 
    [legend_handle,res1, res2] = simu_resi_f1f2(A, B, C, D, B_f, D_f, B_d, D_d, K, P{i}, f1s, f2s, dist, s, x0, x_target);
    title(subplot_titles{i},'Interpreter','latex','FontSize',15);
    [FAR_1, FDR_1] = FARFDR(res1);
    [FAR_2, FDR_2] = FARFDR(res2);
    set(legend_handle, 'Visible', 'off');

        FAR1(1,i) = FAR_1;
        FDR1(1,i) = FDR_1;
        FAR2(1,i) = FAR_2;
        FDR2(1,i) = FDR_2;
end
axis(ha(8), 'off');
set(legend_handle, 'Visible', 'on');
legend_handle.Units = 'normalized';
pos = ha(8).Position;
legend_handle.Position = [pos(1) 0.17 pos(3) 0.2];
legend_handle.FontSize=15;

figure(4)
ha = tight_subplot(2, 4, [0.2 0.05], [0.12 0.07], [0.05 0.05]);
set(4, 'Position', [500, 100, 1200, 450]);

for i = 8:14
    axes(ha(i-7)); 
    [legend_handle,res1, res2] = simu_resi_f1f2(A, B, C, D, B_f, D_f, B_d, D_d, K, P{i}, f1s, f2s, dist, s, x0, x_target);
    title(subplot_titles{i},'Interpreter','latex','FontSize',15);
    [FAR_1, FDR_1] = FARFDR(res1);
    [FAR_2, FDR_2] = FARFDR(res2);
    set(legend_handle, 'Visible', 'off');

        FAR1(2,i-7) = FAR_1;
        FDR1(2,i-7) = FDR_1;
        FAR2(2,i-7) = FAR_2;
        FDR2(2,i-7) = FDR_2;
end
axis(ha(8), 'off');
set(legend_handle, 'Visible', 'on');
legend_handle.Units = 'normalized';
pos = ha(8).Position;
legend_handle.Position = [pos(1) 0.17 pos(3) 0.2];
legend_handle.FontSize=15;


%% Detectability vs Rho
beta_list = [0.01 0.03 0.05 0.1];
omega_list = [10^-2 10^-1.5 10^-1 10^-0.5 10^0 10^0.5 10 10^1.5 100];
color = {'#0072BD','#EDB120','#D95319','#7E2F8E'};
symbol = {'x--','*--','o--','d--'};
theta_list=[theta1,theta2];
for rho = 1:2
P = FD_DRO(A, B, C, D, B_f, D_f, B_d, D_d, s, xi_hat, theta_list(rho), alpha, rho);
detect = caclu_detectablity(V,P,rho);
figure(10+rho*2-1);
plot([0.01 100],[detect detect], 'r--', 'LineWidth', 1);  
hold on;

    for i = 1:size(beta_list,2)
        beta = beta_list(i);
        detect_list = [];
        kappa_list = [];
        for omega = omega_list
            theta = theta_list(rho);            
            [P ,kappa, ~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta, alpha, beta, omega, rho);
            detect_list = [detect_list,caclu_detectablity(V,P,rho)];
            kappa_list = [kappa_list,kappa];
        end
        figure(10+rho*2-1);
        plot(omega_list,detect_list,symbol{i}, 'LineWidth', 1,'color',color{i})
        hold on;
        grid on;
        xlabel('$\omega$', 'Interpreter', 'latex');
        set(gca, 'XScale', 'log')
        if rho == 1
            ylabel('Metric $\rho_1(\cdot)$', 'Interpreter', 'latex');
        else
            ylabel('Metric $\rho_2(\cdot)$', 'Interpreter', 'latex');
        end
        
        figure(10+rho*2);
        plot(omega_list,kappa_list,symbol{i}, 'LineWidth', 1,'color',color{i})
        hold on;
        grid on;
        xlabel('$\omega$', 'Interpreter', 'latex');
        ylabel('$\kappa$', 'Interpreter', 'latex');
        set(gca, 'XScale', 'log')
    end
end


figure(11)
legend({'DR$1$','$\beta=0.99$', '$\beta=0.97$','$\beta=0.95$', ...
    '$\beta=0.90$'}, 'Interpreter', 'latex', 'Location', 'southwest');
set(11, 'Position', [100, 100, 560, 210]);

figure(13)
legend({'DR$2$','$\beta=0.99$', '$\beta=0.97$','$\beta=0.95$', ...
    '$\beta=0.90$'}, 'Interpreter', 'latex', 'Location', 'southwest');
set(13, 'Position', [100, 100, 560, 210]);


figure(12)
    legend({'$\beta=0.99$', '$\beta=0.97$','$\beta=0.95$', ...
    '$\beta=0.90$'}, 'Interpreter', 'latex', 'Location', 'northwest');
    set(12, 'Position', [100, 100, 560, 210]);

figure(14)
    legend({'$\beta=0.99$', '$\beta=0.97$','$\beta=0.95$', ...
    '$\beta=0.90$'}, 'Interpreter', 'latex', 'Location', 'northwest');
    set(14, 'Position', [100, 100, 560, 210]);




