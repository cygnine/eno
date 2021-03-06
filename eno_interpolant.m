function[u] = eno_interpolant(x,y,z,varargin)
% eno_interpolant -- evaluates eno finite-difference interpolant
%
% [u] = eno_interpolant(x,y,z,{k=3})
%
%     Interpolates the data set (x,y) using a piecewise k-th order polynomial
%     using the ENO stencil-choosing rubric, evaluates at the points z. We allow
%     x to be non-equispaced. 
% 
%     The nodal vector x needs to be sorted, but need not contain nodes at the
%     boundaries. If points z lie outside the cells defined by x, extrapolation
%     from the nearest cell is used.

persistent strict_inputs eno_setup newton_eval
if isempty(strict_inputs)
  from labtools import strict_inputs;
  from eno import eno_setup;
  from speclab.newton_polynomials import newton_evaluate as newton_eval
end

opt = strict_inputs({'k'}, {3},[],varargin{:});

% Force column vector
x = x(:);
y = y(:);
zsize = size(z);
z = z(:);

n = length(x);
eno_info = eno_setup(x,y,'k',opt.k);

% Determine indicators for locations of nodes
% Temporarily redefine x as the bin separators to include all real numbers
x= [-inf; x(2:(end-1)); inf];
[temp,bin] = histc(z,x);

u = zeros(size(z));

% For loops are bad, but I can't figure out how to vectorize this
% Compute matrix of locations to interpolate to
for q = 1:n
  flags = bin==q;
  if any(flags)
    u(flags) = newton_eval(eno_info.XInput(q,:),eno_info.dd(:,q),z(flags));
  end
end

u = reshape(u, zsize);
