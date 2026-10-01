function setup_paths()
%SETUP_PATHS  Put lib/ and bundled third-party code on the MATLAB path.
%
% Graph-cuts need max_flow_mex in third_party/matlab_bgl/private/
% (e.g. max_flow_mex.mexw64 on Windows). Do not add private/ to the path.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));
addpath(fullfile(here, 'third_party', 'hernando'));
addpath(fullfile(here, 'third_party', 'matlab_bgl'));

mex_dir = fullfile(here, 'third_party', 'matlab_bgl', 'private');
if isempty(dir(fullfile(mex_dir, 'max_flow_mex.*'))) ...
        || isempty(dir(fullfile(mex_dir, ['max_flow_mex.' mexext])))
    warning('SPARCQ:MissingMaxFlowMex', ...
        ['No max_flow_mex.%s in\n  %s\n' ...
         'Copy your working MEX there (same file you always paste).'], ...
        mexext, mex_dir);
end
end
