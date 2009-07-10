function[indicators] = sobolev_indicators(nc,x);
% [INDICATORS] = SOBOLEV_INDICATORS(NC,X)
%
%     Given Newton polynomial modal coefficients (NC) defining polynomials with
%     basis functions given by the Newton nodal locations (X), computes the
%     highest nontrivial Sobolev norm over every non-overlapping subinterval.
%     I.e. For every column of NC (and also X), computes the Sobolev norms over
%     intervals specified by sequential pairwise subsets of rows of X.  The
%     inputs (DS) are the outputs from taylor_weights_periodic: they are the
%     linear Taylor relations connecting lower-order Taylor expandsions to
%     higher-order expansions


global handles;
newton = handles.speclab.NewtonPolynomials;
mono = handles.speclab.monomials;


k = size(nc,1) - 1; % polynomial order
C = size(nc,2);     % number of `elements'

indicators = zeros([k,C]);
dx = diff(x,[],1);
dx2 = dx.^2;

mc = newton.newton_to_monomial(nc,x);

% This is ridiculously slow and stupid: for each subinterval, for each
% derivative. When I get less lazy: 'l', 'll' are horrible varnames.
for ll = 1:k
  % Select subinterval
  interval = x(ll:(ll+1),:);
  % Reset to original modal coefficients and dx at each pass
  derivative_mc = mc;
  dx_temp = dx(ll,:);

  for l = 1:k
    % Compute next derivative
    derivative_mc = mono.monomial_derivative(derivative_mc);
    % Square it, integrate over interval, add to previous value
    indicators(ll,:) = indicators(ll,:) + dx_temp.*...
        mono.monomial_integrate(mono.monomial_square(derivative_mc),interval);
    % Rescale dx weight to account for next derivative
    dx_temp = dx_temp.*dx2(ll,:);
  end
end
