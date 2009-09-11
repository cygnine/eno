function[stencil,r,indicator] = eno_stencil(x,y,varargin)
% eno_stencil -- returns the stencil necessary for eno operations
%
% [stencil,{r,}] = eno_stencil(x,y,{k:3})
%
%     Given inputs (x,y) that are nodal locations and evaluations, respectively,
%     uses the ENO reconstruction motivation to adaptively choose the least
%     oscillatory stencil for K-th order interpolation. K must be greater than
%     0. If x is of length N, the returned matrix has N stencil indicators: one
%     for each interval of reconstruction between the points. 
%
%     [X(1), X(2)] <---> stencil(1,;)
%     [X(2), X(3)] <---> stencil(2,:)
%                    .
%                    .
%                    .
%     [X(N-1), X(N)] <---> stencil(N-1,:) 
%     
%     The output stencil is the finite-difference stencil used for computing
%     divided differences, and the optional output R is the interval shift
%     relative to `default' stencil interval. 

global handles
cm = handles.common;
fd = handles.FiniteDifference;

opt = cm.InputSchema({'k'}, {3}, [],varargin{:});
k = opt.k;

n = length(x);

% Trival cases
if k==0
  %varargout{1} = zeros([n,1], 'int8');
  %varargout{2} = 0;
  r = zeros([n,1], 'int8');
  indicator = 0;
  stencil = [(1:n).'];
  return 
elseif k==1
  stencil_periodicity = zeros([n,2],'int8');
  stencil_periodicity(n,2) = 1;
  %varargout{1} = stencil_periodicity;
  %varargout{2} = 0;
  r = stencil_periodicity;
  indicator = 0;
  stencil = [(1:n).', mod((2:n+1).'-1,n)+1];
  return 
end

xmin = min(x)-1;  % Completely random values
xmax = max(x)+1;

%%%% Do the same thing as eno_stencil_periodicty
% Must now extend x and y for periodicity: introduce ghost points
x = x(:);
y = y(:);
XGhost = [xmin - (xmax - x(n-k+2:n)); ...
     x; ...
     xmax + x(1:k) - xmin];
YGhost = [y(n-k+2:n); ...
     y; ...
     y(1:k)];

% Chooses stencils in the same fashion as common/FiniteDifference/DifferenceStencil (r=0)
NegativeCount = zeros([n,1],'int32');
PositiveCount = zeros([n,1],'int32');
differences = YGhost;

% Go ahead and pick first point to the right: this is an artificial bias.
RightX = XGhost(2:end);
LeftX = XGhost(1:(end-1));
differences = diff(differences)./(RightX-LeftX);

% These are indices for the differences vector, which has a dynamic size
LeftIndices = int32(((k-1):(k+n-2)).');
RightIndices = int32(((k):(k+n-1)).');

% The strategy here is (surprise) inefficient: computing divided differences
% only to choose stencil.
% TODO: Cherry-pick differences as needed: will make things *much* faster.
for q = 2:k
  % Compute all differences
  RightX = XGhost((q+1):end);
  LeftX = XGhost(1:(n+2*(k-1)-q+1));
  differences = diff(differences)./(RightX-LeftX);

  % Iterate Positive/NegativeCount by doing horrible implicit typecasting
  inds = abs(differences(LeftIndices))<=abs(differences(RightIndices));
  NegativeCount = NegativeCount + int32(inds);
  LeftIndices = LeftIndices - int32(inds);
  RightIndices = RightIndices - int32(inds);
  inds = not(inds);
  PositiveCount = PositiveCount + int32(inds);
end

% Stencil shifts relative to `default' stencil
r = zeros([n,1],'int32');
[PositiveCount, NegativeCount];

r = (PositiveCount-NegativeCount + mod(k,2))/2;
if mod(k,2)==1
  r = (PositiveCount-NegativeCount)/2 + mod(k,2);
else
  r = (PositiveCount-NegativeCount+1)/2;
end

[stencil] = fd.difference_stencil(n,k,'r',r);

%varargout{1} = r;
%varargout{2} = 1;
indicator = 1;
