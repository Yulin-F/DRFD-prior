function [P] = FD_initial(A, B, C, D, B_f, D_f, B_d, D_d, s, f1, xi_hat, theta, alpha)
% Fault detection use initial solution P0

fs = repmat(f1, s + 1, 1);   
[~, V] = matrix_WV(A, B, C, D, B_f, D_f, B_d, D_d, s);

[P_list,~,status] = initial_Iteration(V,theta,alpha,xi_hat,fs) ;
if status ~= 0
    disp('initial no sulotion')
    return
end

P = P_list{end} ;

end 

%% initial_Iteration: Sequential minimization
function [P_list,obj_list,status] = initial_Iteration(V,theta,alpha,xi_hat,f) 
status = 0;
cycle = 100;
P = 0.03 * eye(size(V,1));
eps = 1e-3;
P_list={};
obj_list = [];
for i = 1:cycle
    disp(['initial Iteration ',num2str(i)])
    [tau, sol1] = dro_solve_tau(P,alpha,xi_hat);
    [P, obj, sol2] = initial_solve_P(tau,V,theta,alpha,xi_hat,f);
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

%% initial_solve_P: main problem solve P
function [P_value,obj_value,sol_status] = initial_solve_P(tau,V,theta,alpha,xi_hat,f) 

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

sol = optimize(Constraints, -f'*V'*P*V*f, options);
Objective = -f'*V'*P*V*f;

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


