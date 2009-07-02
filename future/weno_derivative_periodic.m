% MATLAB File : weno_derivative_periodic.m
% [du,{Dks,ds}] = weno_derivative_periodic(x,u,k,interval)
%
% * Creation Date : 2009-06-06
%
% * Last Modified : Fri 12 Jun 2009 03:47:00 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given the nodal locations (x), the function values (u), the ENO
% order of approximation (k), and the interval of approximation (interval),
% calculates the derivative using the WENO approximation assuming periodicity
% over the interval. The nodes x must be ordered. The WENO approximation is
% (2k)-th order accurate in smooth regimes and k-th order `ENO accurate' in
% non-smooth regimes. 
% The optional outputs Dks and ds are, respectively, the k-th order Taylor
% differentiation matrices that are linearly combined to give the WENO
% approximation, and the Taylor-Taylor connection weights connecting the
% aforementioned differentiation matrices to the centered (2k)-th order
% differentiation matrix. These outputs are given as mesh-dependent factors;
% they can be re-used for a fixed mesh x to generate the WENO smoothness
% indicators needed to form a differentiation approximation with a different
% function u on the same mesh.

function[du,varargout] = weno_derivative_periodic(x,u,k,interval)

[Dks,betas,ds] = calculate_weno_weights_periodic(x,u,k,interval);

N = length(u);
du = zeros([N,1]);
for q = 1:(k+1);
  du = du + (spdiags(betas(:,q),0,N,N)*Dks{q})*u;
end

varargout{1} = Dks;
varargout{2} = ds;
