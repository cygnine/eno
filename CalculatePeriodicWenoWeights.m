% MATLAB File : CalculatePeriodicWenoWeights.m
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Sat 06 Jun 2009 03:39:53 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given the grid x and the finite-difference values u, computes the
%   required WENO Taylor weights and Sobolev Indicators. Returns everything
%   separately: the Dks are the individual unweighted lower-order Taylor
%   differentiation matrices, the betas are the already-scaled weights ready for
%   application for these values of u, and the ds are the unscaled Taylor-Taylor
%   connections that can be utilized when computing betas for a different vector
%   u.

function[Dks,betas,ds] = CalculatePeriodicWenoWeights(x,u,k,interval)

global common;
prevpath = addpaths(common.bases.d1.newton.base);

[Dks,ds] = CalculatePeriodicTaylorWeights(x,k,interval);

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
nc = DividedDifference(xinterp,u(inds));

[betas] = SobolevIndicators(nc,xinterp,xinterp([1,k+1],:),ds).';

path(prevpath);

% The betas are applied in the following manner:
%D = zeros(size(Dks{1}));
%for q = 1:(k+1);
%  Dks{q} = spdiags(betas(:,q),0,n,n)*Dks{q};
%  D = D + Dks{q};
%end
% The full differentiation matrix is given by D
