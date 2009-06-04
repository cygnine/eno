% MATLAB File : EnoDerivativePeriodic.m
% [d] = EnoDerivativePeriodic(x,y,k,interval)
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Thu 04 Jun 2009 04:01:37 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given inputs (x,y) that are nodal locations and evaluations,
%   respectively, uses the ENO reconstruction motivation to adaptively choose
%   the least oscillatory stencil for k-th order differentiation. k must be
%   greater than 0. interval is a 2-vector specifying the periodicity of the
%   interval.

function[d] = EnoDerivativePeriodic(x,y,k,interval)

global common;
prevpath = addpaths(common.FiniteDifference);

n = length(x);
xmin = interval(1); xmax = interval(2);

% Must now extend x and y for periodicity: introduce ghost points
x = x(:);
XGhost = [xmin - (xmax - x(n-k+1:n)); ...
     x; ...
     xmax + (x(1:k) - xmin)];
YGhost = [y(n-k+1:n); ...
     y; ...
     y(1:k)];

% Chooses stencils in the same fashion as common/FiniteDifference/DifferenceStencil (r=0)
NegativeCount = zeros([n,1],'int8');
PositiveCount = zeros([n,1],'int8');
differences = YGhost;
dd = zeros([n+2*k,k+1]);
dd(:,1) = YGhost;
% These are indices for the differences vector, which has a dynamic size
LeftIndices = int8((k:(k+n-1)).');
RightIndices = int8(((k+1):(k+n)).');
% The strategy here is (surprise) inefficient: computing divided differences
% only to choose stencil.
% TODO: Cherry-pick differences as needed: will make things *much* faster.
for q = 1:k
  % Compute all differences
  RightX = XGhost((q+1):n+2*k);
  LeftX = XGhost(1:(n+2*k-q));
  differences = diff(differences)./(RightX-LeftX);
  dd(1:(n+2*k-q),q+1) = differences;

  % Iterate Positive/NegativeCount by doing horrible implicit typecasting
  inds = abs(differences(LeftIndices))<=abs(differences(RightIndices));
  NegativeCount = NegativeCount + int8(inds);
  LeftIndices = LeftIndices - int8(inds);
  inds = not(inds);
  PositiveCount = PositiveCount + int8(inds);
  %RightIndices = RightIndices + inds;
end

% Stencil shifts relative to `default' stencil
r = zeros([n,1],'int8');
r = PositiveCount-NegativeCount + mod(k,2);
[stencil,StencilPeriodicity] = DifferenceStencil(n,k,r,true);

% Compute x values
XInput = x(stencil);
inds = StencilPeriodicity==1;
% For indices that wrap down to 1:
XInput(inds) = xmax + (XInput(inds) - xmin);

inds = StencilPeriodicity==-1;
% For indices that wrap up to n:
XInput(inds) = xmin - (xmax - XInput(inds));

% Use stencil to compute interpolants
dd = DividedDifference(XInput.',y(stencil.'));

% Differentiate and evaluate the interpolants
d = NewtonDiffEval(XInput.',dd).';

path(prevpath);
