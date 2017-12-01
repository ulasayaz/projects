
function [out] = generative_model(z, Ws)
  % Simple feed-forward neural network with ReLU
  output = cell(length(Ws),1);
  out = z;
  for i = 1 : length(Ws)
    out = ReLU(Ws{i} * out);
    output{i} = out;
  end
  
end

function output = ReLU(input)
  output = input .* (input > 0);
end

