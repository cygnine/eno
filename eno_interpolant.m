function[u] = eno_interpolant(x,y,z,varargin)
% [U] = ENO_INTERPOLANT(X,Y,Z,{K=3})
%
%     Interpolates the data set (X,Y) using a piecewise K-th order polynomial
%     using the ENO stencil-choosing rubric, evaluates at the points Z. We allow
%     X to be non-equispaced. 
% 
%     The nodal vector X needs to be sorted, but need not contain nodes at the
%     boundaries. 

global handles;
cm = handles.common;
eno = handles.eno;
newton = handles.speclab.NewtonPolynomials;

% Force column vector
x = x(:);
y = y(:);

opt = cm.InputSchema({'k'}, {3},[],varargin{:});
k = opt.k;

% Compute eno stencil
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
  u(flags) = newton.newton_evaluate(z(flags),dd(:,q),XInput(q,:));
end
