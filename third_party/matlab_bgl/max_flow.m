function [flowval cut R F] = max_flow(A,u,v,varargin)
% MAX_FLOW Compute the max flow on A from u to v.
% David Gleich, Stanford University, 2006-2008

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('algname', 'push_relabel','fix_diag',1);
options = merge_options(options, varargin{:});

if options.fix_diag, A = A - diag(diag(A)); end
if check, check_matlab_bgl(A,struct('noneg',1,'nodiag',1)); end

if ~trans, A = A'; end

n = size(A,1);

if nargout == 2
    [flowval cut] = max_flow_mex(A,u,v,lower(options.algname));
elseif nargout >= 3
    [flowval cut ri rj rv] = max_flow_mex(A,u,v,lower(options.algname));
    R = sparse(ri,rj,rv,n,n);
    if ~trans
        R = R';
    end
else
    flowval = max_flow_mex(A,u,v,lower(options.algname));
end

if nargout >= 4
    F = A - R;
end
