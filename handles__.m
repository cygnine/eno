function[hs,pathadditions] = handles__()
% handles__ -- constructs path tree for the eno module
%
% [hs,pathadditions] = handles__()
%
%     Returns directory pointers for common module in HS. PATHADDITIONS is a
%     cell array with a string in each element indicated paths to add to the
%     global path structure. 

% This is by default
hs.base = fileparts(mfilename('fullpath'));

pathadditions = cell(0);
