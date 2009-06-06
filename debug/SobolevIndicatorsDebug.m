% MATLAB File : SobolevIndicatorsDebug.m
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 05 Jun 2009 10:29:57 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging SobolevIndicators

cd ..

%%%%%%%%%%%%%%%% Test polynomials
nc = [[1;1;1;0], [0;2;-1;3]];
x = [[0;1;2;3],[0;1;2;3]];

% poly1: 1 + x + x*(x-1)
%      = x^2 + 1
% dpoly1 = 2*x
% ddpoly1 = 2
% poly2: 2*x - x*(x-1) + 3*x*(x-1)*(x-2)
%      = 3*x^3 - 10*x^2 + 9*x
% dpoly2 = 9*x^2 - 20*x + 9
% ddpoly2 = 18*x - 20
% dddpoly2 = 18

interval = [[0;3],[0;3]];

betas = SobolevIndicators(nc,x,interval);

cd debug;
