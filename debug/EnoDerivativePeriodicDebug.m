% MATLAB File : EnoDerivativePeriodicDebug.m
%
% * Creation Date : 2009-06-04
%
% * Last Modified : Thu 04 Jun 2009 04:01:52 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging EnoDerivativePeriodic

clear
cd ..

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%%%%%%%%%%% Test smooth function
if true
  f = @(x) sin(x);
  df = @(x) cos(x);

  xmin = 0;
  xmax = 2*pi;
  N = 100;
  k = 5;

  x = sort((xmax-xmin)*rand([N,1]) + xmin);

  fx = f(x);
  dfx = df(x);

  d = EnoDerivativePeriodic(x,fx,k,[xmin,xmax]);

  fprintf('Error is %3.3e\n', norm(dfx-d));

end

%%%%%%%%%%%%% Test discontinuous function
if true
  f = @(x) double(x<3);
  df = @(x) 0*x;

  xmin = 0;
  xmax = 2*pi;
  N = 100;
  k = 4;

  x = sort((xmax-xmin)*rand([N,1])+xmin);

  fx = f(x);
  dfx = df(x);
  
  [stencil,r] = EnoStencilPeriodic(x,fx,k,[xmin,xmax]);
  d = EnoDerivativePeriodic(x,fx,k,[xmin,xmax]);

  fprintf('Error is %3.3f\n', norm(dfx-d));
end

cd debug
