% MATLAB File : eno_fd_interpolant_periodic.m
%
% * Creation Date : 2009-06-12
%
% * Last Modified : Fri 12 Jun 2009 03:34:40 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Using the inputted point values (x,y) (i.e. 'finite-difference'),
% reconstructs an interpolating polynomial via the eno construction procedure
% and returns the evaluation of this reconstruction at the point z. The function
% is assumed periodic over the interval `interval'. The optional input k
% determines the order of the interpolant.

function[z] = eno_fd_interpolant_periodic(x,y,z,interval,varargin)

global common
prevpath = addpaths(common.bases.d1.newton.base);

opt = InputSchema({'k'}, {3},[],varargin{:});
k = opt.k

[stencil,stencil_periodicity] = ...
  eno_stencil_interval_periodic(x,y,interval,varargin{:});

xmax = interval(2); xmin = interval(1);
% Compute x values
XInput = x(stencil);
inds = stencil_periodicity==1;
% For indices that wrap down to 1:
XInput(inds) = xmax + (XInput(inds) - xmin);

inds = stencil_periodicity==-1;
% For indices that wrap up to n:
XInput(inds) = xmin - (xmax - XInput(inds));

n = length(x);
% Use stencil to compute interpolants
if k==0
  dd = reshape(y,[1,n]);
else
  dd = divided_difference(XInput.',y(stencil.'));
end

% Determine indicators for locations of nodes
[temp,bin] = histc(z,[x;xmax + x(1) - xmin]);

% For loops are bad, but I can't figure out how to vectorize this
% Compute matrix of locations to interpolate to
for q = 1:n
  flags = bin==q;
  z(flags) = newton_evaluate(z(flags),dd(:,q),XInput(q,:));
end

path(prevpath);
