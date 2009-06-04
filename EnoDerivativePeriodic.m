% MATLAB File : EnoDerivativePeriodic.m
% [d] = EnoDerivativePeriodic(x,y,k)
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 08:39:09 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given inputs (x,y) that are nodal locations and evaluations,
%   respectively, uses the ENO reconstruction motivation to adaptively choose
%   the least oscillatory stencil for k-th order differentiation. k must be
%   greater than 0.

function[d] = EnoDerivativePeriodic(x,y,k)

global common;
prevpath = addpaths(common.FiniteDifference);

n = length(x);

% Chooses stencils in the same fashion as common/FiniteDifference/DifferenceStencil (r=0)
stencils = zeros([n,2],'int8');
NegativeCount = zeros([n,1],'int8');
PositiveCount = zeros([n,1],'int8');
differences = y;

% The strategy here is (surprise) inefficient: computing divided differences
% only to choose stencil.
% Assume we need k nodes to the left and the right
for q = 1:k
  % Compute left difference
  TempStencil = RectifyIndices((1:n)' - NegativeCount - 1);
  TempDifference = 
end

path(prevpath);

end

% Uses periodicity to "right" the indices like 0, n+1
function[indices] =  RectifyIndices(indices,n)
  indices = mod(indices-1,n)+1;
