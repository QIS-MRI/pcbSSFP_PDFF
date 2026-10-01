function old_default = set_matlab_bgl_default(varargin)
% SET_MATLAB_BGL_DEFAULT  David Gleich, Stanford University, 2006-2008

persistent default_options;
if ~isa(default_options,'struct')
    default_options = struct('istrans', 0, 'nocheck', 0, 'full2sparse', 0);
end

if nargin == 0
    old_default = default_options;
else
    old_default = default_options;
    default_options = merge_options(default_options,varargin{:});
end
