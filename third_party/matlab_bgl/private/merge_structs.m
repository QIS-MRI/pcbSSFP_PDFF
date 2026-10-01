function S = merge_structs(A, B)
% MERGE_STRUCTS  David Gleich, 19 October 2005

S = A;
fn = fieldnames(B);
for ii = 1:length(fn)
    if (~isfield(A, fn{ii}))
        S.(fn{ii}) = B.(fn{ii});
    end
end
