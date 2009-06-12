% MATLAB File : calculate_taylor_weights_periodic.m
% [Dks,ds] = calculate_taylor_weights_periodic(x,k,interval)
%
% * Creation Date : 2009-06-04
%
% * Last Modified : Fri 12 Jun 2009 03:23:22 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Calculates the linear WENO weights (not smoothness indicators) to
%   relate (k+1) kth-order Taylor expansions to a (2k)th order Taylor expansion.
%   This is basically a Taylor polynomial calculation for unstructured meshes.
%   Although the weights are indeed calculated, this goes ahead and returns the
%   affinely mapped k-th order differentiation matrices.

function[Dks,ds] = calculate_taylor_weights_periodic(x,k,interval)

global common;
prevpath = addpaths(common.FiniteDifference.base);

% Let's do this the slow and stupid way: linear determination of the weights
% given the (2k)th and k-th order expansion coefficients. 

% The (2k)th coeffs
D2k = derivative_matrix_periodic(x,2*k,interval);
if mod(k,2)==0
  shifts = -(k/2):(k/2);
else
  shifts = -((k-1)/2):((k+1)/2);
end

% The k-th coeffs
Dks = {};
rcount = 1;
for r = shifts
  Dks{rcount} = derivative_matrix_periodic(x,k,interval,r);
  rcount = rcount + 1;
end

n = length(x);
ds = zeros([n,k+1]);

KCount = 0;
NegativeCount = 0;
PositiveCount = 0;
NegativeK = true;

% Do some fancy stuff to determine the weights
for q = 0:(k-1)
  if NegativeK

    temp = 0;
    for r = 1:NegativeCount
      temp = temp - ds(:,r).*...
             full(diag(circshift(Dks{r},[0,k-r])));
    end
    temp = temp + full(diag(circshift(D2k,[0,k-NegativeCount])));
    temp2 = full(diag(circshift(Dks{NegativeCount+1},[0,k-NegativeCount])));
    ds(:,NegativeCount+1) = temp./temp2;
    NegativeCount = NegativeCount + 1;

  else

    temp = 0;
    for r = 1:PositiveCount
      temp = temp-ds(:,r).*...
             full(diag(circshift(Dks{k+2-r},[0,r-k])));
    end
    temp = temp + full(diag(circshift(D2k,[0,PositiveCount-k])));
    temp2 = full(diag(circshift(Dks{k+1-PositiveCount},[0,PositiveCount-k])));
    ds(:,k+1-PositiveCount) = temp./temp2;
    PositiveCount = PositiveCount + 1;

  end
  
  NegativeK = not(NegativeK);
end

if NegativeK
  ds(:,NegativeCount+1) = 1 - sum(ds(:,setdiff(1:(k+1),NegativeCount+1)),2);
else
  ds(:,k+1-PositiveCount) = 1 - sum(ds(:,setdiff(1:(k+1),k+1-PositiveCount)),2);
end

%for q = 1:(k+1)
%  Dks{q} = diag(ds(:,q))*Dks{q};
%end

path(prevpath);
