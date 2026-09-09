function s = merge(s1, s2)
%MERGE Summary of this function goes here
%   Detailed explanation goes here
if ~strcmp(s1.type, s2.type)
    error('charge and current cannot merge together ');
end

s.type = s1.type;

s.pos = [s1.pos; s2.pos];
s.w = [s1.w; s2.w];

if strcmp(s.type, 'charge')
    s.q = [s1.q; s2.q];
elseif strcmp(s.type, 'current')
    s.Idl = [s1.Idl; s2.Idl];
end

end

