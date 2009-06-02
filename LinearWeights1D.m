% Returns the 1D linear weights as defined in [1].
% [1]: Essentially Non-Oscillatory and Weighted Essentially
% Non-Oscillatory Schemes for Hyperbolic Conservation Laws, C.-W. Shu, 1997
function[crj] = LinearWeights1D(x,k);

% k \in {1,2,3,4,...}
% Assumes a stencil [x_{-(k-1)}, ..., x_{-2}, x_{-1}, x_0,x_1,x_2,..., x_{k-1}] 
% The point at which these constants compute derivatives is x_0

% Allocation
crj = zeros([k+1,k]);

% For each value of r:
for r = -1:(k-1)
  for j = 0:(k-1)
    % c(r+2,j+1)
    for m = (j+1):k
      den = 1
      for l = 0:k
        if l ~= m
          den = den*x(
