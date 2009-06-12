% MATLAB File : eno_fd_interpolant.m
% [uz] = EnoFDInterpolant(x,ux,z,[k=3])
%
% * Creation Date : 2009-06-11
%
% * Last Modified : Fri 12 Jun 2009 03:34:28 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Evaluates the ENO interpolant at z given the finite-difference
%   point-values (x,ux) of a function. Can handle any order k>0, and unstructured
%   meshes as well (though in this case the reconstruction is not conservative.)

function[uz] = EnoFDInterpolant(x,ux,z,varargin)

disp("you have to implement this");
uz = [];
