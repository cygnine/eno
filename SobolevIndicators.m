% MATLAB File : SobolevIndicators.m
% [betas] = SobolevIndicators(nc,x,interval)
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 05 Jun 2009 10:32:00 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given Newton polynomial modal coefficients (nc) defining
%   polynomials with basis functions given by the Newton nodal locations (x),
%   computes the highest nontrivial Sobolev norm over the specified interval
%   (interval).

function[betas] = SobolevIndicators(nc,x,interval);

global common;
prevpath = addpaths(common.bases.d1.newton.base, ...
                    common.bases.d1.monomial.base);

k = size(nc,1) - 1; % polynomial order
C = size(nc,2);     % number of `elements'

betas = zeros([1,C]);
dx = diff(interval,[],1);
dx2 = dx.^2;

mc = NewtonToMonomial(nc,x);

for l = 1:k;
  mc = MonomialDerivative(mc);
  betas = betas + dx.*MonomialIntegrate(MonomialSquare(mc),interval);
  dx = dx.*dx2;
end

path(prevpath);
