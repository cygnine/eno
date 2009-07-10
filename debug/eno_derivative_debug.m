% Debugging eno_stencil

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

  d = eno_derivative(x,fx,x,'k',k);

  fprintf('Error is %3.3e\n', norm(dfx-d));

end

%%%%%%%%%%%%% Test discontinuous function
if true
  f = @(x) double(x<1) + double(x>3);
  df = @(x) 0*x;

  fx = f(x);
  dfx = df(x);
  
  [stencil,r] = eno_stencil(x,fx,'k',k);
  d = eno_derivative(x,fx,x,'k',k);

  fprintf('Error is %3.3f\n', norm(dfx-d));
end

cd debug
