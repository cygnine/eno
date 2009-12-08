function[p] = eno_reconstruction(x,y,varargin)
% eno_reconstruction -- constructs an ENO PiecewisePolynomial
%
% [p] = eno_reconstruction(x,y,{k=3,alpha=0,beta=0,interval=[]})
%
%     Interpolates the data set (x,y) using a piecewise k-th order polynomial
%     using the ENO stencil-choosing rubric. Takes the d'th derivative. See
%     eno_interpolant and eno_derivative.
%
%     The output p is a PiecewisePolynomial object with cell boundaries x. The
%     polynomial is the ENO reconstruction, which may be evaluated anywhere.
%
%     If the optional input interval is specified, it is assumed that the
%     reconstruction should be periodic. 
%
%     The basis_representation for p is the Jacobi polynomial basis of class
%     (alpha,beta).

persistent strict_inputs repnodes eno enop gq
if isempty(strict_inputs)
  from labtools import strict_inputs
  from piecewise_interpolation.grid_tools import replicate_local_nodes as repnodes
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as gq
  from eno import eno_interpolant as eno
  from eno import eno_interpolant_periodic as enop
end

inputs = {'k', 'alpha', 'beta', 'interval'};
defaults = {3, 0, 0, []};
opt = strict_inputs(inputs, defaults, [], varargin{:});

x = x(:);
y = y(:);

[r,w] = gq(opt.k+1,'alpha',opt.alpha, 'beta', opt.beta);

if isempty(opt.interval) % No periodicity
  z = repnodes(r, x);
  fz = eno(x,y,z,'k',opt.k);

  temp.cell_boundaries = x;
else  % Assuming length(interval)==2, it's periodic
  assert(opt.interval(1)<=x(1), 'Input interval must enclose data x');
  assert(opt.interval(2)>=x(end), 'Input interval must enclose data x');

  if (x(1) - opt.interval(1) + opt.interval(2) - x(end))==0  
    % Then ignore y(end)
    z = repnodes(r, x);
    fz = enop(x(1:(end-1)), y(1:(end-1)), z, opt.interval, 'k', opt.k);
    temp.cell_boundaries = x;
  else
    if (x(1) - opt.interval(1))==0
      % Then append interval(2) to end of cell markers
      tempx = [x; opt.interval(2)];
    elseif (opt.interval(2) - x(end))==0
      % Then append interval(1) to beginning of cell markers
      tempx = [opt.interval(1); x];
    else
      % Append both sides of interval to cell markers
      tempx = [opt.interval(1); x; opt.interval(2)];
    end
    z = repnodes(r, tempx);
    fz = enop(x, y, z, opt.interval, 'k', opt.k);
    temp.cell_boundaries = tempx;
  end
end

temp.N = opt.k+1;
temp.x = z;
temp.y = fz;
temp.jacobi_alpha = opt.alpha;
temp.jacobi_beta = opt.beta;

p = PiecewisePolynomial(temp);
