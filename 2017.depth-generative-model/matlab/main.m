close all; clear; clc;

%% example parameters
k = 5;
weight_type = 'gaussian'; % 'gaussian','[0,1]','{0,1}','{-1,1}'
measurement_type = 'subsample'; % 'gaussian','subsample'
% weight_dim = [50, 100, 200, 400];
weight_dim = [50, 100, 200];
% m = 10 * k;
m = weight_dim(end) / 2;

%% iterate
count_success = 0; % count how many times we find the global minimum
iterations = 1000;
epsilon = 1e-2;
total_time = 0;
results_holder = cell(iterations, 1);
for i = 1:iterations
  if mod(i, 100) == 0
    fprintf('==> iteration %05i: \n', i)
  end
  results = toy_problem(k, m, weight_dim, weight_type, measurement_type, false);
  results_holder{i} = results;
  total_time = total_time+results.time;
  if norm(results.z_hat - results.z_gt) < epsilon
    count_success = count_success + 1;
  else
    fprintf('iter %05i: failed to find global minimum. \n\tz_gt =[%s] \n\tz_hat_1=[%s] \n\tz_hat_2=[%s]\n', i, ...
      arr2str(results.z_gt), arr2str(results.z_hat_1), arr2str(results.z_hat_2));
  end
end
fprintf('success rate=%g%% \n', 100.0 * count_success / iterations)
fprintf('mean(time)=%g \n', total_time / iterations)











