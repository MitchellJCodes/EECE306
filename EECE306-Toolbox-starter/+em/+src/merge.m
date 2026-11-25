function s = merge(s1, s2, varargin)
% Merge compatible discretized sources.
% Requires two sources and accepts additional sources.

% Check whether the sources are the same type
if ~strcmp(s1.type, s2.type)
    error('Charge and current cannot merge together.');
end

s.type = s1.type;

% Merge positions and weights
s.pos = [s1.pos; s2.pos];
s.w = [s1.w(:); s2.w(:)];

% Merge source elements
if strcmp(s.type, 'charge')

    s.q = [s1.q(:); s2.q(:)];
    s.Idl = [];

elseif strcmp(s.type, 'current')

    s.Idl = [s1.Idl; s2.Idl];
    s.q = [];

else

    error('Invalid source type.');

end

s.label = '';

% Merge any additional sources
for k = 1:numel(varargin)
    s = em.src.merge(s, varargin{k});
end

end
