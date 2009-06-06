% MATLAB File : WenoDifferentiationDebug.m
%
% * Creation Date : 2009-06-06
%
% * Last Modified : Sat 06 Jun 2009 03:46:25 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging/testing for CalculatePeriodicWenoWeights, primarily
%   aimed at differentiation

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
cd ..

%%%%%%%%%%% Try arctan function
if true
  xmin = 0;
  xmax = 1;
  interval = [xmin,xmax];
  k = 4;
  C = 15;
  f = @(x) atan(C*(x-(xmax-xmin)/2));
  df = @(x) C./(1+(C*(x-(xmax-xmin)/2)).^2);
  N = 100;

  x = sort((xmax-xmin)*rand([N,1]) + xmin);
  %x = linspace(xmin,xmax,N+1).';
  %x = x(1:N);
  u = f(x);
  du = WenoDerivativePeriodic(x,u,k,interval);
  
  EnoDu = EnoDerivativePeriodic(x,u,k,interval);

  fprintf('Eno error is %3.3e\nWeno error is %3.3e\n', norm(EnoDu-df(x)),...
          norm(du-df(x)));
  semilogy(x,abs(df(x) - EnoDu), 'b', x, abs(du - df(x)), 'r')
  
end

cd debug
