% Tests the eno stencil evaluation procedure

global handles;

N = 100;
x = 2*rand([N-1,1]) - 1;
x = sort([x;-1]);
interval = [-1,1];
k = 3;
C = 25;

f = @(x) atan(C*x) - x*atan(C);
df = @(x) C./(1+(C*x).^2);

fx = f(x);

cd ..
[stencil,stencil_periodicity,reference_offset] = ...
         eno_stencil_periodic(x,fx,interval,'k',k);

z = linspace(-1,1,1e3).';

u = handles.eno.eno_interpolant_periodic(x,fx,z,interval,'k',k);
du = handles.eno.eno_derivative_periodic(x,fx,z,interval,'k',k) + atan(C);

cd debug
