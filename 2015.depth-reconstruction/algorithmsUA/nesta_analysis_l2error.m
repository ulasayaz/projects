function [x,timeNesterov,iter] = nesta_analysis_l2error(D,b,A,samples,epsilon)

%% initialize

Nm2 = size(D,1);
N = size(D,2);
nonSamples = setdiff([1:N],samples);

Nabla_f_mu = zeros(N,1);
u_mu = zeros(Nm2,1); 
Nabla_f_mu_cum = zeros(N,1);

iter = 1000;

mu = 0.01;

for k = 1:iter
    
    alpha = 1/(2*(k+1));
    tau = 2/(k+3);
    
    Dx = D * x;
    ind_dx_small = find(abs(Dx) < mu);
    u_mu(ind_dx_small) = 1/mu * Dx(ind_dx_small);
    ind_dx_large = setdiff(allEntries,ind_dx_small);
    u_mu(ind_dx_large) = sign(Dx(ind_dx_large));
    Nabla_f_mu = D' * u_mu;
    
    Nabla_f_mu_cum = Nabla_f_mu_cum + alpha * Nabla_f_mu;
    
    
    
end


