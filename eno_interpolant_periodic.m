function[u] = eno_interpolant_periodic(x,y,z,interval,varargin)
% [U] = ENO_INTERPOLANT_PERIODIC(X,Y,Z,INTERVAL,{K=3})
%
%     Interpolates the data set (X,Y) using a piecewise K-th order polynomial
%     using the ENO stencil-choosing rubric, evaluates at the points Z. We allow
%     X to be non-equispaced.  The default is K=3.
% 
%     This is a periodic extension, so the 2-vector INTERVAL specifying the
%     domain is required. The nodal vector X needs to be sorted, but need not
%     contain nodes at the boundaries. 

global handles;
cm = handles.common;
eno = handles.eno;
newton = handles.speclab.NewtonPolynomials;

% Force column vector
x = x(:);
y = y(:);

opt = cm.InputSchema({'k'}, {3},[],varargin{:});
k = opt.k;

[stencil,stencil_periodicity] = eno.eno_stencil_periodic(x,y,interval,'k',k);

xmax = interval(2); xmin = interval(1);
% Compute x values
XInput = x(stencil);
inds = stencil_periodicity==1;
% For indices that wrap down to 1:
XInput(inds) = xmax + (XInput(inds) - xmin);

inds = stencil_periodicity==-1;
% For indices that wrap up to n:
XInput(inds) = xmin - (xmax - XInput(inds));

n = length(x);
% Use stencil to compute interpolants
if k==0
  dd = reshape(y,[1,n]);
else
  dd = newton.divided_difference(XInput.',y(stencil.'));
end

% Determine indicators for locations of nodes
[temp,bin] = histc(z,[x;xmax + x(1) - xmin]);

u = zeros(size(z));

% For loops are bad, but I can't figure out how to vectorize this
% Compute matrix of locations to interpolate to
for q = 1:n
  flags = bin==q;
  %u(flags) = newton.newton_evaluate(z(flags),dd(:,q),XInput(q,:));
  u(flags) = newton.newton_evaluate(XInput(q,:),dd(:,q),z(flags));
end
