function [P] = FD_DRO(A, B, C, D, B_f, D_f, B_d, D_d, s, xi_hat, theta, alpha, rho)
% Fault detection use DRFD

[~, V] = matrix_WV(A, B, C, D, B_f, D_f, B_d, D_d, s);

[P_list,~,status] = dro_Iteration(V,theta,alpha,xi_hat,rho) ;
if status ~= 0
    disp('DRO no sulotion')
    return
end

% figure(99);
% plot(obj_list, 'b-o', 'MarkerSize', 5, 'LineWidth', 1.5);
% if rho == 1
%     xlabel('$\rho_1(P) + \gamma \eta$', 'Interpreter', 'latex','FontSize', 15);
% else
%     xlabel('$\rho_2(P) + \gamma \eta$', 'Interpreter', 'latex','FontSize', 15);
% end
% ylabel('Cross-validation FAR', 'Interpreter', 'latex','FontSize', 15);
% grid on

P = P_list{end} ;

end

%% dro_Iteration: Sequential minimization
function [P_list,obj_list,status] = dro_Iteration(V,theta,alpha,xi_hat,rho) 
status = 0;
cycle = 100;
P = 0.01 * eye(size(V,1));
eps = 1e-5;
P_list={};
obj_list = [];
for i = 1:cycle
    disp(['DRO Iteration ',num2str(i)])
    [tau, sol1] = dro_solve_tau(P,alpha,xi_hat);
    [P, obj, sol2] = dro_solve_P(tau,V,theta,alpha,xi_hat,rho)   ;
    P_list{i} = P;
    obj_list(end +1) = obj;
    if i > 1
        if abs(obj_list(end) -  obj_list(end - 1)) <= eps
            break
        end
        if sol1 ~= 0 || sol2 ~= 0
            disp("iteration no solution");
            status = 1;
            return
        end       
    end
   
end

end

%% dro_solve_P: main problem solve P
function [P_value,obj_value,sol_status] = dro_solve_P(tau,V,theta,alpha,xi_hat,rho) 

N = size(xi_hat,2);

P = sdpvar(size(V,1),size(V,1),'symmetric');
y = sdpvar(N, 1);
t = sdpvar(N, 1);
q = sdpvar(N, 1);
lambda = sdpvar(1, 1);


Constraints = [
    theta + (1/N) * sum(y) - lambda * alpha <= 0;
    y >= 0;
    y >= lambda * ones(N,1) - t;
    lambda >= 0;
    P>=0;
];

for i = 1:N
    Constraints = [Constraints; 
                   [eye(size(P))-tau(i)*P, -xi_hat(:,i);
                   -xi_hat(:,i)', xi_hat(:,i)'*xi_hat(:,i)-q(i)+tau(i)] >= 0;]; 

    Constraints = [Constraints; 
                   [q(i), t(i);
                    t(i), 1] >= 0;]; 

end


options = sdpsettings('verbose', 1, 'solver', 'mosek','verbose',0, 'mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS', 1000, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 3600);

if rho == 1 
sol = optimize(Constraints, -trace(V'*P*V), options);
Objective = -trace(V'*P*V);
elseif rho == 2
    [U1, S, ~] = svd(V, 'econ');
sol = optimize(Constraints, -logdet(S'*U1'*P*U1*S), options);
Objective = -logdet(S'*U1'*P*U1*S);
else
    disp("rho can be only chosen as 1 or 2");
end
P_value = value(P);
obj_value = value(Objective);
sol_status = sol.problem;
end

%% dro_solve_tau: subproblem J_1 solve tau
function [tau_value,sol_status] = dro_solve_tau(P,alpha,xi_hat)

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
    Objective == (1/N) * sum(y) - lambda * alpha;];

for i = 1:N
    Constraints = [Constraints; 
                   [eye(size(P))-tau(i)*P, -xi_hat(:,i);
                   -xi_hat(:,i)', xi_hat(:,i)'*xi_hat(:,i)-q(i)+tau(i)] >= 0;]; 

    Constraints = [Constraints; 
                   [q(i), t(i);
                    t(i), 1] >= 0;]; 
end

options = sdpsettings('verbose', 0, 'solver', 'mosek','verbose',0, 'mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS', 1000, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 3600);
sol = optimize(Constraints, Objective, options);

tau_value = value(tau);
sol_status = sol.problem;
end