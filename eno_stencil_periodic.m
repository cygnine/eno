% MATLAB File : eno_stencil_periodic.m
% function[stencil,{r}] = eno_stencil_periodic(x,y,k,interval)
%
% * Creation Date : 2009-06-04
%
% * Last Modified : Fri 12 Jun 2009 03:36:21 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given inputs (x,y) that are nodal locations and evaluations,
%   respectively, uses the ENO reconstruction motivation to adaptively choose
%   the least oscillatory stencil for k-th order differentiation. k must be
%   greater than 0. interval is a 2-vector specifying the periodicity of the
%   interval.
%   
%   The output stencil is the finite-difference stencil used for computing
%   divided differences, and the optional output r is the interval shift
%   relative to `default' differentiation interval. There is no assumption made
%   about where the nodes lie; the reconstruction as if all the action is
%   happening at exactly the grid points. 

function[stencil,varargout] = eno_stencil_periodic(x,y,k,interval)

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
NegativeCount = zeros([n,1],'int32');
PositiveCount = zeros([n,1],'int32');
differences = YGhost;
% These are indices for the differences vector, which has a dynamic size
LeftIndices = int32((k:(k+n-1)).');
RightIndices = int32(((k+1):(k+n)).');
% The strategy here is (surprise) inefficient: computing divided differences
% only to choose stencil.
% TODO: Cherry-pick differences as needed: will make things *much* faster.
for q = 1:k
  % Compute all differences
  RightX = XGhost((q+1):n+2*k);
  LeftX = XGhost(1:(n+2*k-q));
  differences = diff(differences)./(RightX-LeftX);

  % Iterate Positive/NegativeCount by doing horrible implicit typecasting
  inds = abs(differences(LeftIndices))<=abs(differences(RightIndices));
  NegativeCount = NegativeCount + int32(inds);
  LeftIndices = LeftIndices - int32(inds);
  %%%%%%%%%%%%% Is this line necessary?
  RightIndices = RightIndices - int32(inds);
  %%%%%%%%%%%%%
  inds = not(inds);
  PositiveCount = PositiveCount + int32(inds);
end

% Stencil shifts relative to `default' stencil
r = zeros([n,1],'int32');
r = (PositiveCount-NegativeCount + mod(k,2))/2;
[stencil,stencil_periodicity] = difference_stencil(n,k,r,true);

varargout{1} = r;
