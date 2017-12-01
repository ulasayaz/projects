close all; clear; clc;

% Experiment settings
settings.k = 10;
settings.weight_type = 'gaussian'; % 'gaussian','[0,1]','{0,1}','{-1,1}'
settings.measurement_type = 'gaussian'; % 'gaussian','subsample'
settings.iterations = 500;
settings.m = 40;
n1_array = settings.k:5:80;
n2_array = settings.m:5:120;

% Compute the success rate
success_array = zeros(length(n1_array), length(n2_array));
for i = 1 : length(n1_array)
  n1 = n1_array(i);
  parfor j = 1 : length(n2_array)
    n2 = n2_array(j);
    success_array(i,j) = success_vs_n1_n2(n1, n2, settings);
    fprintf('n1=%03d, n2=%03d, success rate=%.2f \n', n1, n2, success_array(i,j))
  end
end

% Visualization of the n1-n2 phase transition plot
X = repmat(n1_array', 1, length(n2_array));
Y = repmat(n2_array, length(n1_array), 1);
surf(X, Y, success_array);
xlabel('n1')
ylabel('n2')
zlabel('success rate [%]')
title(sprintf('k=%d, m=%d, sampling=%s, iter=%d', ...
  settings.k, settings.m, settings.measurement_type, settings.iterations))
hold on;
contour(X, Y, success_array, 'LevelStep', 5)

% Function for computing the success rate of a (n1-n2) pair
function success_rate = success_vs_n1_n2(n1, n2, settings)
count_success = 0; % count how many times we find the global minimum
weight_dim = [n1, n2];
epsilon = 1e-2;
total_time = 0;
for i = 1 : settings.iterations
  results = toy_problem(settings.k, settings.m, weight_dim, ...
    settings.weight_type, settings.measurement_type, false);
  total_time = total_time+results.time;
  if norm(results.z_hat - results.z_gt) < epsilon
    count_success = count_success + 1;
  end
end
success_rate = 100.0 * count_success / settings.iterations;
end
