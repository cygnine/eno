function[eno] = init__()
% init__ -- Initialization file for eno package
%
% [nodes] = init__()

eno = recurse_files(pwd);
eno.debug = matlab_import('debug');
eno.sceno = matlab_import('sceno');
