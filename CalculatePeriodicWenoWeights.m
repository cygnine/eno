% MATLAB File : CalculatePeriodicWenoWeights.m
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 05 Jun 2009 10:58:56 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given the grid x and the finite-difference values u, computes the
%   required WENO Taylor weights and Sobolev Indicators. Returns these
%   separately; they can be combined as needed, and the Sobolev indicators can
%   be replaced when u changes but the mesh remains the same.

function[Dks,betas] = CalculatePeriodicWenoWeights(u,x,k);
