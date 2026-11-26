function s = merge(varargin)
%MERGE Concatenate compatible discretized sources.
%   S = EM.SRC.MERGE(S1,S2,...) concatenates the element lists of two or
%   more sources of the same type ('charge' or 'current').
%
%   Each source must contain positions and weights, plus q for charge
%   sources or Idl for current sources. Mixing source types is an error.

if nargin < 2
    error('em.src.merge:TooFewInputs', ...
        'At least two source arguments are required.');
end

s = varargin{1};
if ~isstruct(s) || ~isfield(s,'type') || ...
        ~isfield(s,'pos') || ~isfield(s,'w')
    error('em.src.merge:InvalidSource', ...
        'Each source must be a struct with type, pos, and w fields.');
end
if ~any(strcmp(s.type, {'charge','current'}))
    error('em.src.merge:InvalidType', ...
        'Source type must be ''charge'' or ''current''.');
end

for k = 2:nargin
    sk = varargin{k};
    if ~isstruct(sk) || ~isfield(sk,'type') || ...
            ~isfield(sk,'pos') || ~isfield(sk,'w')
        error('em.src.merge:InvalidSource', ...
            'Each source must be a struct with type, pos, and w fields.');
    end
    if ~strcmp(s.type, sk.type)
        error('em.src.merge:MixedTypes', ...
            'Cannot merge a %s source with a %s source.', s.type, sk.type);
    end
    if ~any(strcmp(sk.type, {'charge','current'}))
        error('em.src.merge:InvalidType', ...
            'Source type must be ''charge'' or ''current''.');
    end
    s.pos = [s.pos; sk.pos];
    s.w = [s.w(:); sk.w(:)];
    if strcmp(s.type,'charge')
        if ~isfield(s,'q') || ~isfield(sk,'q')
            error('em.src.merge:MissingChargeData', ...
                'Charge sources must contain q.');
        end
        s.q = [s.q(:); sk.q(:)];
    else
        if ~isfield(s,'Idl') || ~isfield(sk,'Idl')
            error('em.src.merge:MissingCurrentData', ...
                'Current sources must contain Idl.');
        end
        s.Idl = [s.Idl; sk.Idl];
    end
end

if ~isfield(s,'label')
    s.label = '';
end
end
