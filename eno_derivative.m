function[u] = eno_derivative(x,y,z,varargin)
% eno_derivative -- computes an eno finite-difference derivative
%
% [u] = eno_derivative(x,y,z,{k=3, d=1})
%
%     Computes the d'th derivative of the ENO interpolant at the nodal locations Z.
%     See ENO_INTERPOLANT. Interpolates the data set (X,Y) using a
%     piecewise K-th order polynomial using the ENO stencil-choosing rubric,
%     takes the derivative, and evaluates at the points Z. We allow X to be
%     non-equispaced.  
% 
%     The nodal vector X needs to be sorted, but need not contain nodes at the
%     boundaries. 

global packages;
cm = packages.labtools;
eno = packages.eno;
newton = packages.speclab.newton_polynomials;

% Force column vector
x = x(:);
y = y(:);
zsize = size(z);
z = z(:);

opt = cm.input_schema({'k', 'd'}, {3, 1},[],varargin{:});
k = opt.k;

stencil = eno.eno_stencil(x,y,'k',k);

% Compute x values
XInput = x(stencil);

n = length(x);
% Use stencil to compute interpolants
if k==0
  dd = reshape(y,[1,n]);
else
  dd = newton.divided_difference(XInput.',y(stencil.'));
end

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
    %u(flags) = newton.newton_derivative_evaluate(XInput(q,:).',dd(:,q),z(flags));
    u(flags) = newton.newton_derivative_evaluate(XInput(q,:).',dd(:,q),'z',z(flags), 'd', opt.d);
  end
end

u = reshape(u,zsize);
