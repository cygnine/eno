function[u] = speno_interpolation(x,y,z,varargin)
% speno_interpolation -- evaluate speno finite-difference interpolant
%
% [u] = speno_interpolation(x,y,z,{k=3})
%
%     Interpolates the data set (x,y) using a piecewise k-th order polynomial
%     using the SPENO stencil-choosing rubric, evaluates at the points z. We allow
%     x to be non-equispaced. 
% 
%     The nodal vector x needs to be sorted, but need not contain nodes at the
%     boundaries. If points z lie outside the cells defined by x, extrapolation
%     from the nearest cell is used.

persistent newton_eval strict_inputs eno_setup bisection gq
if isempty(newton_eval)
  from speclab.newton_polynomials import newton_evaluate as newton_eval
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as gq
  from labtools import strict_inputs;
  from labtools.rootfind import bisection
  from eno import eno_setup
end

opt = strict_inputs({'k'}, {3},[],varargin{:});

% Force column vector
x = x(:);
y = y(:);
zsize = size(z);
z = z(:);

n = length(x);
eno_info = eno_setup(x,y,'k',opt.k);

% Now we need to determine if there are any stencil choices that warrant
% intra-cell shock placement. Note that this is impossible to do unless there
% are k buffer cells on either side from the boundary.
shock_flags = false([max([1,n-1-2*opt.k]),1]);

% Determine extremal value of the shift r:
rmax = ceil(opt.k/2);
rmin = rmax - (opt.k-1);

% Step across stencil: find cells where neighboring reconstructions avoid the
% cell. These cells probably contain shocks.
cellmin = opt.k+1;
cellmax = (n-1) - opt.k;
cell_select = cellmin:cellmax;
left_leaning = (eno_info.r(cell_select-1) == rmin);
right_leaning = (eno_info.r(cell_select+1) == rmax);

shock_flags = left_leaning & right_leaning;

%%%%%%% 
% But how to deal with situations when the shock cell falsely sets off a
% neighboring cell as harboring a shock?
%
% Well, I guess we'll just assume this never happens
%%%%%%%

% Define a function on the shock cell to be the right interpolant minus the left
% interpolant. If the function has a root in the shock cell, we'll place that
% shock. Otherwise, we make a discontinuity explicit.
N_shocks = sum(shock_flags);  % Number of possible shocks
shock_id = find(shock_flags);
shock_id = shock_id + cellmin - 1;  % remember that shock_flags elides some cells
local_nodes = gq(opt.k+1);
offset = 0;   % a counter for how many shocks we place
for q = 1:N_shocks;
  shock_cell = shock_id(q);
  a = x(shock_cell-offset);
  b = x(shock_cell-offset+1);
  shift = (a+b)/2; scale = (b-a)/2;

  % First define a polynomial that is the difference of the interpolants
  cell_nodes = local_nodes*scale + shift;

  f = @(x) newton_eval(eno_info.XInput(shock_cell-offset+1,:), eno_info.dd(:,shock_cell-offset+1), x) ...
         - newton_eval(eno_info.XInput(shock_cell-offset-1,:), eno_info.dd(:,shock_cell-offset-1), x);

  [intersection, bisection_flag] = bisection(a, b, f);

  if bisection_flag==4  % f(a)*f(b) > 0, so there wasn't a single simple root
    % In 1D, I don't think there's a need to do anything; just keep the
    % interpolant as it is
  else  % I guess we found an intersection; basically remove the shock cell
    x(shock_cell-offset) = intersection;
    x(shock_cell-offset+1) = [];
    eno_info.XInput(shock_cell+offset,:) = [];
    eno_info.dd(:,shock_cell+offset) = [];

    offset = offset + 1;
  end
end

%%%%%%%%%%%%% 
%  Carry on with usual eno interpolation routine
%%%%%%%%%%%%% 

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
