% Tests the eno stencil evaluation procedure

global handles;

N = 100;
x = 2*rand([N,1]) - 1;
x = sort(x);
interval = [-1,1];
k = 3;

f = @(x) abs(x)<0.5;

fx = f(x);

cd ..
[stencil,stencil_periodicity,reference_offset] = ...
         eno_fd_stencil_periodic(x,fx,interval,k);

z = linspace(-1,1,1e3);

u = handles.eno.eno_fd_interpolant_periodic(x,fx,z,interval,k);

cd debug
