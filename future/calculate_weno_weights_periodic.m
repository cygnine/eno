function[Dks,betas,ds] = calculate_weno_weights_periodic(x,u,k,interval)
% [DKS,BETAS,DS] = CALCULATE_WENO_WEIGHTS_PERIODIC(X,U,K,INTERVAL)

%   Given the grid X and the finite-difference values U, computes the
%   required WENO Taylor weights and Sobolev Indicators. Returns everything
%   separately: the DKS are the individual unweighted lower-order Taylor
%   differentiation matrices, the BETAS are the already-scaled weights ready for
%   application for these values of U, and the DS are the unscaled Taylor-Taylor
%   connections that can be utilized when computing betas for a different vector
%   U.

global common;
prevpath = addpaths(common.bases.d1.newton.base);

[Dks,ds] = calculate_taylor_weights_periodic(x,k,interval);

n = length(x);
inds = zeros([k+1,n]);
PositiveFlag = false([k+1,n]);
NegativeFlag = false([k+1,n]);
for q = 1:(k+1)
  inds(q,:) = (1:n)+(q-1);
end

PositiveFlag = inds>n;
NegativeFlag = inds<1;
inds = mod(inds-1,n)+1;

% Calculate modal coefficient for Newton polynomial basis
xinterp = x(inds);
xinterp(PositiveFlag) = xinterp(PositiveFlag) + (interval(2) - interval(1));
xinterp(NegativeFlag) = xinterp(NegativeFlag) - (interval(2) - interval(1));
nc = divided_difference(xinterp,u(inds));

[betas] = sobolev_indicators(nc,xinterp,xinterp([1,k+1],:),ds).';

path(prevpath);

% The betas are applied in the following manner:
%D = zeros(size(Dks{1}));
%for q = 1:(k+1);
%  Dks{q} = spdiags(betas(:,q),0,n,n)*Dks{q};
%  D = D + Dks{q};
%end
% The full differentiation matrix is given by D
