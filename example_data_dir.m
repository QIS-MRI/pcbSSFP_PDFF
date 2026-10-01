function d = example_data_dir()
%EXAMPLE_DATA_DIR  example data/  (created if missing)

d = fullfile(fileparts(mfilename('fullpath')), 'example data');
if ~isfolder(d)
    mkdir(d);
end
end
