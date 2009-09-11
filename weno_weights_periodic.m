function[Dks,beta_matrix,ds] = weno_weights_periodic(x,u,k,interval)
% weno_weights_periodic -- periodic weno weight calculation
%
% [dks,betas,ds] = weno_weights_periodic(x,u,k,interval)
%
%     Given the grid X and the finite-difference values U, computes the required
%     WENO Taylor weights and Sobolev Indicators. Returns everything separately:
%     the DKS are the individual unweighted lower-order Taylor differentiation
%     matrices, the BETAS are the already-scaled weights ready for application
%     for these values of U.

global handles;
eno = handles.eno;
newton = handles.speclab.NewtonPolynomials;

[Dks,ds] = eno.taylor_weights_periodic(x,k,interval);

% Hack-y way of determining what the stencil is for each interpolant. I'm sure
% there's a better way.
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

% Calculate modal coefficients for Newton polynomial basis
xinterp = x(inds);
xinterp(PositiveFlag) = xinterp(PositiveFlag) + (interval(2) - interval(1));
xinterp(NegativeFlag) = xinterp(NegativeFlag) - (interval(2) - interval(1));
% Compute modal coefficients from (x,y) points
nc = newton.divided_difference(xinterp,u(inds));

% Use local modal coefficients to compute sobolev norms
[sobolev_norms] = eno.sobolev_indicators(nc,xinterp);

epsilon = 1e-6;  % tolerance
betas = 1./(epsilon+sobolev_norms).^2;

% Time for fun index games: match-up globally indexed Sobolev norms with locally
% indexed stencil references
beta_matrix = zeros([k,n]);
for q = 1:k;
  beta_matrix(q,:) = circshift(betas(k+1-q,:),[0,k-q]).*ds(:,q).';
end
beta_matrix = beta_matrix*spdiags(1./sum(beta_matrix,1).',0,n,n);

% The betas are applied in the following manner for differentiation:
%  
%  D = zeros(size(Dks{1}));
%  for q = 1:(k+1);
%    Dks{q} = spdiags(betas(:,q),0,n,n)*Dks{q};
%    D = D + Dks{q};
%  end
%
% The full differentiation matrix is given by D
