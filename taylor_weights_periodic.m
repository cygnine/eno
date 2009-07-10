function[Dks,ds] = calculate_taylor_weights_periodic(x,k,interval)
%
% [DKS,DS] = CALCULATE_TAYLOR_WEIGHTS_PERIODIC(X,K,INTERVAL)
%
%     Calculates the linear WENO weights (not smoothness indicators) to relate K
%     Kth-order Taylor expansions to a (2K-1)th order Taylor expansion.  This is
%     basically a Taylor polynomial calculation for unstructured meshes.
%     Although the weights are indeed calculated, this goes ahead and returns
%     the affinely mapped K-th order differentiation matrices.
% 
%     The details are as follows: you give me an unstructured mesh of size N. I
%     return to you an (N+1) x K set of coefficients. For row m of that matrix,
%     the K coefficients correspond to the linear combination of K-th order
%     stencils needed to return a (2K-1)-th order interpolation stencil. The
%     column location reflects the offset of the K-th order stencil. Row m
%     corresponds to the interval between X(m) and X(m+1), where X(m+1) == X(1)
%     if m+1 == N+1. 

global handles;
fd = handles.FiniteDifference;

% Let's do this the slow and stupid way: linear determination of the weights
% given the (2k-1)th and k-th order expansion coefficients. We'll use
% rudimentary linear algebra to determine the weights by matching the
% coefficients for the derivative stencils evaluated at the left-hand point of
% the interval.

% The (2k)th coeffs, (FD diff stencil shifted one to the right)
D2k = fd.derivative_matrix_periodic(x,2*k-1,interval,ones(size(x)));
if mod(k,2)==1
  shifts = (-((k-1)/2):((k-1)/2)) + 1;
else
  shifts = -(k/2-1):(k/2);
end

% The k-th coeffs
Dks = {};
rcount = 1;
for r = shifts
  Dks{rcount} = fd.derivative_matrix_periodic(x,k,interval,r*ones(size(x)));
  rcount = rcount + 1;
end

n = length(x);
ds = zeros([n,k]);

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
             full(diag(circshift(Dks{r},[0,k-1-NegativeCount])));
    end
    temp = temp + full(diag(circshift(D2k,[0,k-1-NegativeCount])));
    temp2 = full(diag(circshift(Dks{NegativeCount+1},[0,k-1-NegativeCount])));
    ds(:,NegativeCount+1) = temp./temp2;
    NegativeCount = NegativeCount + 1;

  else

    temp = 0;
    for r = 1:PositiveCount
      temp = temp-ds(:,k-r+1).*...
             full(diag(circshift(Dks{k+1-r},[0,PositiveCount-k])));
    end
    temp = temp + full(diag(circshift(D2k,[0,PositiveCount-k])));
    temp2 = full(diag(circshift(Dks{k-PositiveCount},[0,PositiveCount-k])));
    ds(:,k-PositiveCount) = temp./temp2;
    PositiveCount = PositiveCount + 1;

  end
  
  NegativeK = not(NegativeK);
end

if NegativeK
  ds(:,NegativeCount) = 1 - sum(ds(:,setdiff(1:k,NegativeCount)),2);
else
  ds(:,k-PositiveCount) = 1 - sum(ds(:,setdiff(1:k,k-PositiveCount)),2);
end
