function[u] = sceno_interpolation(x,y,z,varargin)
% sceno_interpolation -- evaluate sceno finite-difference interpolant
%
% [u] = sceno_interpolation(x,y,z,{k=3})
%
%     Interpolates the data set (x,y) using a piecewise k-th order polynomial
%     using the SCENO stencil-choosing rubric, evaluates at the points z. We allow
%     x to be non-equispaced. 
% 
%     The nodal vector x needs to be sorted, but need not contain nodes at the
%     boundaries. If points z lie outside the cells defined by x, extrapolation
%     from the nearest cell is used.

%newton_eval = packages.speclab.newton_polynomials.newton_evaluate.handle;
newton_eval = from_as('speclab.newton_polynomials', 'newton_evaluate');
from labtools input_schema;
from eno eno_setup

opt = input_schema({'k'}, {3},[],varargin{:});

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
%%%%%%%
