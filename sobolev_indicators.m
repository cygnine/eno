function[BetaMatrix] = sobolev_indicators(nc,x,interval,ds);
% [BETAMATRIX] = SOBOLEV_INDICATORS(NC,X,INTERVAL,DS)
%
%     Given Newton polynomial modal coefficients (nc) defining polynomials with
%     basis functions given by the Newton nodal locations (x), computes the
%     highest nontrivial Sobolev norm over the specified interval (interval).
%     The inputs (ds) are the outputs from CalculatePeriodicTaylorWeights: they
%     are the linear Taylor relations connecting lower-order Taylor expandsions
%     to higher-order expansions

epsilon = 1e-6;  % tolerance

global common;
prevpath = addpaths(common.bases.d1.newton.base, ...
                    common.bases.d1.monomial.base);

k = size(nc,1) - 1; % polynomial order
C = size(nc,2);     % number of `elements'

betas = zeros([1,C]);
dx = diff(interval,[],1);
dx2 = dx.^2;

mc = newton_to_monomial(nc,x);

for l = 1:k;
  mc = monomial_derivative(mc);
  betas = betas + dx.*monomial_integrate(monomial_square(mc),interval);
  dx = dx.*dx2;
end

betas = 1./(epsilon+betas).^2;

BetaMatrix = zeros([k+1,C]);
for q = 1:(k+1);
  BetaMatrix(q,:) = circshift(betas,[0,k-(q-1)]).*ds(:,q).';
end
BetaMatrix = BetaMatrix*spdiags(1./sum(BetaMatrix,1).',0,C,C);

path(prevpath);
