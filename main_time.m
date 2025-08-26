clear;clc;close all;
rng(1113);
omega = [0.01 0.1 1 ];
s = 4;  
times = 200000;
alpha = 0.01;
beta = [0.05 0.03 0.01];
sample = 200;
theta1 = 0.008;
theta2 = 0.009;

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
xi_hat = W * [matrix_laplace(size(W,2), sample, 0, var), matrix_laplace(size(W,2), 5000, 0, var)];
N_sample = [10 20 40 80 160 320 640 1280];
theta = [theta1,theta2];
%% FD design
i =1;
for N = N_sample 
   for rho= 1:2
        tic;
        [P1 ,~, ~] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat(:,1:N), theta(rho), alpha, beta(1), omega(3), rho);
        time = toc;
        time_list(rho,i) = time;
    end
    i=i+1;
end

% Time fig 
figure(1)
plot(N_sample,time_list(1,:),'*-','LineWidth', 1);
hold on;
plot(N_sample,time_list(2,:),'o-','LineWidth', 1);

ylabel('Solution time (s)', 'Interpreter', 'latex','FontSize', 15);
xlabel('$N$', 'Interpreter', 'latex','FontSize', 15);
legend({'NDR$1$-C','NDR$2$-C'}, 'Interpreter', 'latex','FontSize', 15,'Location','northwest');
grid on
% set(gca,'XScale','log');
set(gca,'YScale','log');
set(1, 'Position', [100, 100, 550, 250]);