function h = sparcq_waitbar(h, frac, msg)
%SPARCQ_WAITBAR  Create, update, or close the SPARCQ waitbar.
%
%   h = sparcq_waitbar([], 0, 'Starting...')
%   sparcq_waitbar(h, 0.5, 'Halfway')
%   sparcq_waitbar(h, 'close')

if nargin >= 2 && (ischar(frac) || isstring(frac)) && strcmp(frac, 'close')
    if ~isempty(h) && ishghandle(h)
        delete(h);
    end
    h = [];
    return
end

frac = max(0, min(1, double(frac)));
if nargin < 3
    msg = '';
end

if isempty(h) || ~ishghandle(h)
    if isempty(msg)
        msg = 'SPARCQ';
    end
    h = waitbar(frac, msg, 'Name', 'SPARCQ');
    return
end

if isempty(msg)
    waitbar(frac, h);
else
    waitbar(frac, h, msg);
end
end
