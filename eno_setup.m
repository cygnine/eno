function[eno_info] = eno_setup(x,y,varargin)
% eno_setup -- ENO setup calculations
%
% eno_info = eno_divided_differences(x,y,{k=3, interval=false, periodic=false})
%
%     Computes the 0-th through k-th Newton divided differences for the data
%     (x,y), where the stencil is chosen according to the ENO rubric. If
%     `periodic' is true, `interval', must be a 2-element vector specifying the
%     interval of periodicity. Also returns the eno stencil and the periodicity
%     indicators (if appropriate). 

global packages;
inputs = {'k', 'interval', 'periodic'};
defaults = {3, false, false};
opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});
divided_difference = packages.speclab.newton_polynomials.divided_difference.handle;
eno = packages.eno.eno_stencil.handle;
enop = packages.eno.eno_stencil_periodic.handle;

% Force column vectors:
x = x(:);
y = y(:);

if opt.periodic
  if length(opt.interval)<2
    error('If the approximation is periodic, you must specify a bounding interval');
  end
end

if opt.periodic
  [stencil, stencil_periodicity] = enop(x,y,opt.interval, 'k', opt.k);
  interval = opt.interval;
  xmax = interval(2); xmin = interval(1);
  % Compute x values
  XInput = x(stencil);
  inds = stencil_periodicity==1;
  % For indices that wrap down to 1:
  XInput(inds) = xmax + (XInput(inds) - xmin);

  inds = stencil_periodicity==-1;
  % For indices that wrap up to n:
  XInput(inds) = xmin - (xmax - XInput(inds));
  eno_info.stencil_periodicity = stencil_periodicity;
else
  stencil = eno(x,y,'k',opt.k);
  stencil = stencil(1:end-1,:);

  XInput = x(stencil);
end

n = length(x);
% Use stencil to compute interpolants
if opt.k==0
  eno_info.dd = reshape(y,[1,n]);
else
  eno_info.dd = divided_difference(XInput.',y(stencil.'));
end

eno_info.stencil = stencil;
eno_info.XInput = XInput;
