function [obj, grad, hessian] = empirical_risk_obj(x,x0, W, d, A, ts)


xii = x;
for ii = 1:d
    xii = W{ii}*xii;
    Wpx{ii} = W{ii}' * diag(xii >= 0);
    xii = max(0,xii - ts(ii));
end

% paul's old code for computing backprop
% for ii = 1:d   
%     xii = max(0,W{ii}*xii - ts(ii));
%     Wpx{ii} = W{ii}(xii >= 0, :);
% end

xoii = x0;
for ii = 1:d
    xoii = max(0,W{ii}*xoii - ts(ii));
    Wpxo{ii} = W{ii}(xoii >= 0, :);
end

% Compute Objective
lambda = 0.0;
obj = 1/2 * 1/size(A, 1) * norm(A * (xii - xoii))^2 + lambda * norm(x)^2;

% Compute Gradient & Hessian
grad= A'*A*(xii - xoii);
hessian = A';
for ii = d:-1:1
    grad = Wpx{ii}*grad;
    hessian = Wpx{ii}*hessian;
end

grad = grad / size(A,1) + 2 * lambda * x;
hessian = hessian*hessian';

