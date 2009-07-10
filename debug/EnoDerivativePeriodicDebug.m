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
xmin = 0;
xmax = 2*pi;
N = 100;
k = 5;

x = sort((xmax-xmin)*rand([N,1])+xmin);

%%%%%%%%%%%% Test smooth function
if true
  f = @(x) sin(x);
  df = @(x) cos(x);

  fx = f(x);
  dfx = df(x);

  d = eno_derivative_periodic(x,fx,x,[xmin,xmax],'k',k);

  fprintf('Error is %3.3e\n', norm(dfx-d));

end

%%%%%%%%%%%%% Test discontinuous function
if true
  f = @(x) double(x<1) + double(x>3);
  df = @(x) 0*x;

  fx = f(x);
  dfx = df(x);
  
  [stencil,r] = eno_stencil_periodic(x,fx,[xmin,xmax],'k',k);
  d = eno_derivative_periodic(x,fx,x,[xmin,xmax],'k',k);

  fprintf('Error is %3.3f\n', norm(dfx-d));
end

cd debug
