function options = merge_options(default_options,varargin)
% MERGE_OPTIONS  David F. Gleich, Stanford University, 2008

if ~isempty(varargin) && mod(length(varargin),2) == 0
    options = merge_structs(struct(varargin{:}),default_options);
elseif length(varargin)==1 && isstruct(varargin{1})
    options = merge_structs(varargin{1},default_options);
elseif ~isempty(varargin)
    error('matlag_bgl:optionsParsing',...
        'There were an odd number of key-value pairs of options specified');
else
    options = default_options;
end
