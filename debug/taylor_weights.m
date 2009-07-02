% Script for debugging taylor weights

global handles;
fd = handles.FiniteDifference;
eno = handles.eno;

ks = 2:7;
N = 100;
interval=[-1,1];
x = linspace(interval(1),interval(2),N+1);
x = x(1:N).';

for k=ks
  D2k = fd.derivative_matrix_periodic(x,2*k-1,interval,ones(size(x)));
  [Dks,ds] = eno.calculate_taylor_weights_periodic(x,k,interval);

  error = full(D2k);
  for q = 1:k
    error = error - full(spdiags(ds(:,q),0,N,N)*Dks{q});
  end

  fprintf('Norm for k=%1f is %1.6e\n', k, norm(error));
end
