% MATLAB File : CalculatePeriodicWenoWeights.m
% [ds] = CalculatePeriodicWenoWeights(x,k,interval)
%
% * Creation Date : 2009-06-04
%
% * Last Modified : Thu 04 Jun 2009 07:22:06 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Calculates the linear WENO weights (not smoothness indicators) to
%   relate (k+1) kth-order Taylor expansions to a (2k)th order Taylor expansion.
%   This is basically a Taylor polynomial calculation for unstructured meshes.

function[ds] = CalculatePeriodicWenoWeights(x,k,interval)

global common;
prevpath = addpaths(common.FiniteDifference.base);

D2k = DerivativeMatrixPeriodic(x,2*k,interval);
if mod(k,2)==0
  shifts = -(k/2):(k/2);
else
  shifts = -((k-1)/2):((k+1)/2);
end

Dks = {};
rcount = 1;
for r = shifts
  Dks{rcount} = DerivativeMatrixPeriodic(x,k,interval,r);
  rcount = rcount + 1;
end

n = length(x);
ds = zeros([n,k+1]);



path(prevpath);
