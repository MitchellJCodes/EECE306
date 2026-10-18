function s = merge(s1, s2)
%MERGE Concatenation Function
%   S = EM.SRC.MERGE(s1, s2) returns a struct that joins the element lists 
%   of two point sources together.
%   
%   Example:
%   s1 = em.src.pointCharge(1e-9, [-0.1 0 0]);
%   s2 = em.src.pointCharge(1e-9, [ 0.1 0 0]);
%   em.src.merge(s1, s2); % returns a struct with fields:
%   type = 'charge'
%   pos = [-0.1000, 0, 0; 0.1000, 0, 0]
%   w = [1; 1]
%   q = [1.0000e-09; 1.0000e-09]
%
%   See also POINTCHARGE.

if nargin < 2
    error('At least two source arguments are required.');
end

if ~strcmp(s1.type, s2.type)
    error('ERROR: Source types do not match.');
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

