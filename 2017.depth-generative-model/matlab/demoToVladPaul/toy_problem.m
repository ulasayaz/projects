%TOY_PROBLEM( k, m, n1, n2, weight_type, doplot ) a toy problem for depth image reconstruction
%   k: the dimension of the latent input vector
%   m: the dimension of the measurement vector
%   weight_dim: the dimensions of the neural network layer
%   weight_type: distribution of neural network weights, choose from 'gaussian','[0,1]','{0,1}','{-1,1}'
%   measurement_type: type of sampling matrix, choose from 'gaussian','subsample'
%   doplot: do plot or not, only valid when k=2

  k = 2;
  m = k*20;
  weight_dim = [30, 100];
 
  doplot = true; % plots the cost function. set 'true' only if k=2
  
%% Ground Truth
z_gt = randn(k, 1);

%% Generative model
Ws = cell(1,length(weight_dim));
Ws{1} = create_weights(k, weight_dim(1), 'gaussian');
for i = 2 : length(weight_dim)
  Ws{i} = create_weights(weight_dim(i-1), weight_dim(i), 'gaussian');
end
% W1 = create_weights(k, n1, weight_type);
% W2 = create_weights(n1, n2, weight_type);
% W3 = create_weights(n2, n3, weight_type);
% G = @(z)generative_model_simple(z, W1, W2, W3);
G = @(z)generative_model(z, Ws);
x_gt = G(z_gt);

%% Sampling
num_out = weight_dim(end);
A = create_weights(num_out, m, 'gaussian'); 
y = A * x_gt;

%% Optimization
% We are using a two-stage algorithm here
% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
options = optimoptions('fminunc','Display', 'off', 'UseParallel', true);
C = @(z)(norm(A * G(z) - y));
z0 = randn(k, 1);
t = cputime;
z_hat_1 = fminunc(C, z0, options);
z_hat_2 = fminunc(C, -z_hat_1, options);
% z_hat_2 = z_hat_1;
t = cputime - t;
if C(z_hat_1)<C(z_hat_2)
  z_hat = z_hat_1;
else
  z_hat = z_hat_2;
end
x_hat = G(z_hat);


  %% Evaluation
  diff_z = norm(z_hat-z_gt);
  diff_x = norm(x_hat-x_gt);
  diff_C = C(z_hat)-C(z_gt);
  fprintf('|z_hat - z_gt|=%g \n', diff_z)
  fprintf('|x_hat - x_gt|=%g \n', diff_x)
  fprintf('|C_hat - C_gt|=%g \n', diff_C)
  fprintf('z_hat/z_gt=%.3f, z_gt/z_hat=%.3f \n', norm(z_hat)/norm(z_gt), norm(z_gt)/norm(z_hat))

  %% Plot cost function landscape
if doplot
  if (k==2) % && (diff_z > 1e-2)
    boundary = -7:0.02:7;
    [XX,YY] = meshgrid(boundary,boundary);
    XX = XX + z_gt(1);
    YY = YY + z_gt(2);
    P = [XX(:),YY(:)];
    Cs = zeros(size(P,1), 1);
    for i = 1 : size(P,1)
      Cs(i) = C(P(i,:)');
    end
    Cs = reshape(Cs, size(XX));
    
    figure(1);
    cla;
    % surf(XX, YY, Cs);
    mesh(XX, YY, Cs);
    hold on
    contour(XX, YY, Cs, 'LevelStep', max(Cs(:))/100 )
    colormap('jet')
    xlabel('z(1)')
    ylabel('z(2)')
    zlabel('cost function $$||y - AG(z)||$$', 'Interpreter', 'Latex')
    title({'Cost Function Landscape', ...
      sprintf('$$|z_{hat} - z_{gt}|=%g$$', diff_z), ...
      sprintf('$$|G(z_{hat}) - G(z_{gt})|=%g$$', diff_x), ...
      sprintf('$$|C_{hat} - C_{gt}|=%g$$', diff_C)}, ...
      'Interpreter', 'Latex')
    hold on;
    z_offset=max(Cs(:))/100;
    p2=scatter3(z_gt(1), z_gt(2), C(z_gt)+z_offset, 'og', ...
      'MarkerFaceColor', 'g');
    p3=scatter3(z_hat_1(1), z_hat_1(2), C(z_hat_1)+z_offset, 'hm', ...
      'MarkerFaceColor', 'm');
    p4=scatter3(z_hat_2(1), z_hat_2(2), C(z_hat_2)+z_offset, 'hk', ...
      'MarkerFaceColor', 'k');
    
    % contourf(XX, YY, Cs);
    leg = legend([p2, p3, p4], 'ground truth', '$\hat{z_1}$', '$\hat{z_2}$')
    set(leg,'Interpreter','latex');
  else
    disp('perfect reconstruction, so no plotting. program ends here.')
  end
end

function W = create_weights(num_in, num_out, type)
  if strcmp(type,'[0,1]')
    W = rand(num_out, num_in);
  elseif strcmp(type,'gaussian')
    variance = 1 / num_in;
    std = sqrt(variance);
    W = std * randn(num_out, num_in);
  elseif strcmp(type,'{0,1}')
    W = rand(num_out, num_in) > 0.5;
  elseif strcmp(type,'{-1,1}')
    W = 2 * (rand(num_out, num_in) > 0.5) - 1;
  else
    error('invalid type')
  end
end


