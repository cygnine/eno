function[c] = coefficient_representation(x,y,varargin)
% coefficient_representation -- ENO-based functional construction
%
% [c] = coefficient_representation(x,y,{interval=false, k=3, alpha=0, beta=0})
%
%     Uses the ENO stencil-choosing rubric on the point values (x,y) to
%     interpolate a k-th order polynomial on each cell [x(i), x(i+1)].
%
%     If the 2-element vector `interval' is specified, the function is assumed
%     to be periodic on the specified interval. The output c is then a (k+1) x
%     length(x) matrix, where c(:,i) are the (k+1) L^2-normalized
%     Jacobi(alpha,beta) modal coefficients representing the polynomial on cell
%     [x(i), x(i+1)]. The coefficients c(:,end) are the coefficients on the
%     remainder of the periodic interval. Each set of coefficients corresponds
%     to local coordinates on [-1,1].
%
%     If `interval' is not specified, the function is not assumed to have
%     periodicity and c is a (k+1) x (length(x) - 1) matrix.

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
eno_setup = packages.eno.eno_setup.handle;
newton_eval = packages.speclab.newton_polynomials.newton_evaluate.handle;

inputs = {'interval', 'k', 'alpha', 'beta'};
defaults = {false, 3, 0, 0};
opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});

k = opt.k;
x = x(:);
y = y(:);

[r,w] = jac.quad.gauss_quadrature(k+1,'alpha',opt.alpha,'beta',opt.beta);
polys = jac.eval.eval_jacobi_poly(r,0:k,'alpha',opt.alpha,'beta',opt.beta);

cell_scale = diff(x)/2;
cell_shift = (x(2:end) + x(1:(end-1)))/2;
if length(opt.interval)<2
  eno_info = eno_setup(x,y,'k', opt.k);
else
  eno_info = eno_setup(x,y,'k', opt.k, 'periodic', true, 'interval', opt.interval);
  cell_scale(end+1) = ((opt.interval(2) - x(end)) + (x(1) - opt.interval(1)))/2;
  cell_shift(end+1) = x(end) + 1/2*((opt.interval(2) - x(end)) + (x(1) - opt.interval(1)));
end

% Global vertices
N_cell = size(eno_info.stencil,1);
x = repmat(r, [1, N_cell])*spdiags(cell_scale,0,N_cell,N_cell) + ...
    repmat(cell_shift.',[k+1, 1]);

fx = zeros(size(x));

% Evaluate the interpolant at these data points:
for q = 1:N_cell
  fx(:,q) = newton_eval(eno_info.XInput(q,:), eno_info.dd(:,q), x(:,q));
end

% Use the quadrature rule to figure out the modal coefficients
c = polys'*spdiags(w,0,k+1,k+1)*fx;
