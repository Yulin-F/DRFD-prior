function obj = caclu_detectablity(V, P ,rho)

if rho == 1 
obj = trace(V'*P*V);
elseif rho == 2
    [U1, S, ~] = svd(V, 'econ');
obj = log(det(S'*U1'*P*U1*S));
else
    disp("rho can be only chosen as 1 or 2");
end

end
