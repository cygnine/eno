function[u] = eno_interpolant_periodic(x,y,z,interval,varargin)
% eno_interpolant_periodic -- interpolates periodic pointwise data using eno
%
% [u] = eno_interpolant_periodic(x,y,z,interval,{k=3})
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
newton_eval = handles.speclab.newton_polynomials.newton_evaluate.handle;
eno_setup = eno.eno_setup.handle;

% Force column vector
x = x(:);
y = y(:);
zsize = size(z);
z = z(:);

opt = cm.input_schema({'k'}, {3},[],varargin{:});

%[stencil,stencil_periodicity] = eno.eno_stencil_periodic(x,y,interval,'k',k);

xmax = interval(2); xmin = interval(1);

n = length(x);
eno_info = eno_setup(x,y,'periodic', true, 'interval', interval, 'k', opt.k);

% Determine indicators for locations of nodes
[temp,bin] = histc(z,[x;xmax + x(1) - xmin]);

% Take care of periodic nodes < x(1)
flags = (bin==0);
bin(flags) = n;
flags2 = flags & z<=x(1);
z(flags2) = xmax + (z(flags2) - xmin);

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
