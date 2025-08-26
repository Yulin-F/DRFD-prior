function [P ,kappa, obj_list] = FD_DRO_prior(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta, alpha, beta, omega, rho)
% Fault detection use DRFD-prior

K = place(A,B,[0.5 0.3 -0.3]);
[~, V] = matrix_WV(A, B, C, D, B_f, D_f, B_d, D_d, s);
fs = repmat(f1, s + 1, 1);

[P0] = FD_initial(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta, alpha);
[P_list,kappa_list,obj_list,status] = pri_Iteration(V,fs,omega,theta,alpha,beta,xi_hat,P0,rho) ;
if status ~= 0
    disp('RS no sulotion')
    return
end

figure(99);
plot(-obj_list, '-o', 'MarkerSize', 5, 'LineWidth', 1);
if rho == 1
    ylabel('$\rho_1(P) + \gamma \eta$', 'Interpreter', 'latex','FontSize', 15);
else
    ylabel('$\rho_2(P) + \gamma \eta$', 'Interpreter', 'latex','FontSize', 15);
end
xlabel('Iteration', 'Interpreter', 'latex','FontSize', 15);
grid on
set(99, 'Position', [100, 100, 560, 210]); 

P= P_list{end};
kappa = kappa_list(end);

end

%% pri_Iteration: Sequential minimization
function [P_list,kappa_list,obj_list,status] = pri_Iteration(V,f,omega,theta,alpha,beta,xi_hat,P,rho) 
status = 0; 
cycle = 500;
kappa = 1e-5;
eps = 1e-5;
P_list={};
P_list{1} = P;
kappa_list = [kappa];
obj_list = [];
for i = 1:cycle
    disp(['RS Iteration ',num2str(i)])
    [pi, sol1] = pri_solve_pi(kappa,P,xi_hat,V,f);
    [tau, sol2] = pri_solve_tau(P,alpha,xi_hat);
    [ P,kappa,obj, sol3] = pri_solve_P_kappa(tau,pi,V,f,omega,theta,alpha,beta,xi_hat,rho) ;
    P_list{i} = P;
    kappa_list(end +1) = kappa;
    obj_list(end +1) = obj;
    if i > 1
        if abs(obj_list(end) -  obj_list(end - 1)) <= eps
            break
        end
        if sol1 ~= 0 || sol2 ~= 0 || sol3 ~= 0
            disp("iteration no solution");
            status = 1;
            return
        end 
    end
end
end

%% pri_solve_P_kappa: main problem solve P kappa
function [P_value,kappa_value,obj_value, sol_status] = pri_solve_P_kappa(tau,pi,V,f,omega,theta,alpha,beta,xi_hat,rho) 

N = size(xi_hat,2);
Vf = V * f;

%Var in main Problem
P = sdpvar(size(V,1),size(V,1),'symmetric');

kappa = sdpvar(1,1);
%Var in Problem J_1
y = sdpvar(N, 1);
t = sdpvar(N, 1);
q = sdpvar(N, 1);
lambda = sdpvar(1, 1);
%Var in Problem J_2
v = sdpvar(N, 1);
u = sdpvar(N, 1);
p = sdpvar(N, 1);

Constraints = [
    theta + (1/N) * sum(y) - lambda * alpha <= 0;
    (1/N) * sum(v)-beta * kappa<= 0;
    kappa >= 0;
    kappa <= 1e5;
    y >= 0;%Problem J_1
    y >= lambda * ones(N,1) - t;
    lambda >= 0;
    q >= zeros(N,1);
    t >= zeros(N,1);
    v >= 0; %Problem J_2
    v >= kappa * ones(N,1) - u;
    P >= 0;
    p >= zeros(N,1);
    u >= zeros(N,1);
];

for i = 1:N
    Constraints = [Constraints; %Problem J_1
                   [eye(size(P))-tau(i)*P, -xi_hat(:,i);
                   -xi_hat(:,i)', xi_hat(:,i)'*xi_hat(:,i)-q(i)+tau(i)] >= 0;]; 

    Constraints = [Constraints; 
                   [q(i), t(i);
                    t(i), 1] >= 0;]; 

    Constraints = [Constraints; %Problem J_2
                   [eye(size(P))+pi(i)*P, -xi_hat(:,i)+pi(i)*P'*Vf;
                   -xi_hat(:,i)'+pi(i)*Vf'*P, xi_hat(:,i)'*xi_hat(:,i)-p(i)+(Vf'*P*Vf-1)*pi(i)] >= 0;]; 

    
    Constraints = [Constraints; 
                   [p(i), u(i);
                    u(i), 1] >= 0;]; 

end



options = sdpsettings('verbose', 0, 'solver', 'mosek', 'mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS', 10000, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 3600);
if rho == 1 
    sol = optimize(Constraints, -trace(V'*P*V) - omega * kappa, options);
%     Objective = -trace(V'*P*V);
    Objective = -trace(V'*P*V) - omega * kappa;
elseif rho == 2
    [U1, S, U2] = svd(V, 'econ');
    sol = optimize(Constraints, -log(geomean(S'*U1'*P*U1*S)) - omega * kappa, options);   
%     Objective = -logdet(S'*U1'*P*U1*S);
    Objective = -log(geomean(S'*U1'*P*U1*S)) - omega * kappa;
else
    disp("rho can be only chosen as 1 or 2");
end

P_value = value(P);
kappa_value = value(kappa);
obj_value = value(Objective);
sol_status = sol.problem;

end

%% pri_solve_pi: subproblem J_2 solve pi
function [pi_value, sol_status] = pri_solve_pi(kappa,P,xi_hat,V,f)

N = size(xi_hat,2);
Vf = V * f;

v = sdpvar(N, 1);
u = sdpvar(N, 1);
p = sdpvar(N, 1);
pi = sdpvar(N, 1);
Objective = sdpvar(1, 1);

Constraints = [v >= zeros(N,1);
    v >= kappa * ones(N,1) - u;
    pi >= zeros(N,1);
    Objective == (1/N) * sum(v)  ;
    p >= zeros(N,1);
    u >= zeros(N,1);];

for i = 1:N

    Constraints = [Constraints; 
                   [eye(size(P))+pi(i)*P, -xi_hat(:,i)+pi(i)*P'*Vf;
                   -xi_hat(:,i)'+pi(i)*Vf'*P, xi_hat(:,i)'*xi_hat(:,i)-p(i)+(Vf'*P*Vf-1)*pi(i)] >= 0;]; 

    
    Constraints = [Constraints; 
                   [p(i), u(i);
                    u(i), 1] >= 0;]; 
end

options = sdpsettings('verbose', 1, 'solver', 'mosek','verbose',0, 'mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS', 1000, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 3600);


sol = optimize(Constraints, Objective, options);
pi_value = value(pi);
sol_status = sol.problem;

end

%% pri_solve_tau: subproblem J_1 solve tau
function [tau_value, sol_status] = pri_solve_tau(P,alpha,xi_hat)

N = size(xi_hat,2);

y = sdpvar(N, 1);
t = sdpvar(N, 1);
q = sdpvar(N, 1);
tau = sdpvar(N, 1);
lambda = sdpvar(1, 1);
Objective = sdpvar(1, 1);

Constraints = [y >= 0;
    y >= lambda * ones(N,1) - t;
    lambda >= zeros(N,1);
    tau >= zeros(N,1);
    Objective == (1/N) * sum(y) - lambda * alpha;
    q >= zeros(N,1);
    t >= zeros(N,1);];

for i = 1:N
    Constraints = [Constraints; 
                   [eye(size(P))-tau(i)*P, -xi_hat(:,i);
                   -xi_hat(:,i)', xi_hat(:,i)'*xi_hat(:,i)-q(i)+tau(i)] >= 0;]; 

    Constraints = [Constraints; 
                   [q(i), t(i);
                    t(i), 1] >= 0;]; 
end

options = sdpsettings('verbose', 1, 'solver', 'mosek','verbose',0, 'mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS', 1000, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 3600);

sol = optimize(Constraints, Objective, options);
tau_value = value(tau);
sol_status = sol.problem;

end

