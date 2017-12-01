close all; clear; clc;

settings.k = 10;
settings.weight_type = 'gaussian'; % 'gaussian','[0,1]','{0,1}','{-1,1}'
settings.measurement_type = 'subsample'; % 'gaussian','subsample'
% settings.weight_dim = [50, 100, 200, 400];
settings.weight_dim = [50, 100, 200];
% m = 10 * k;
% m = weight_dim(end) / 2;
m_array = [1, settings.k:5:100];
success_array = zeros(1, length(m_array));
iterations = 100;
for i = 1 : length(m_array)
  m = m_array(i);
  success_array(i) = success_vs_measurements(m, iterations, settings);
  fprintf('m=%03d, success rate=%.2f \n', m, success_array(i))
end
plot(m_array, success_array, '-');
xlabel('number of measurements')
ylabel('success rate [%]')
title(sprintf('k=%d, sampling=%s', settings.k, settings.measurement_type))

function success_rate = success_vs_measurements(m, iterations, settings)
count_success = 0; % count how many times we find the global minimum

epsilon = 1e-2;
total_time = 0;
for i = 1:iterations
  results = toy_problem(settings.k, m, settings.weight_dim, ...
    settings.weight_type, settings.measurement_type, false);
  total_time = total_time+results.time;
  if norm(results.z_hat - results.z_gt) < epsilon
    count_success = count_success + 1;
  end
end
success_rate = 100.0 * count_success / iterations;
end










