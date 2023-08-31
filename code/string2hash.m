%STRING2HASH Convert a string to a 64 char hex hash string (256 bit hash)
%
%   Copied from Oliver Woodfords answer in this forum post (5 May 2016):
%   https://se.mathworks.com/matlabcentral/answers/45323-how-to-calculate-hash-sum-of-a-string-using-java#answer_220897
%
function hash = string2hash(string)
persistent md
if isempty(md)
    md = java.security.MessageDigest.getInstance('SHA-256');
end
hash = sprintf('%2.2x', typecast(md.digest(uint8(string)), 'uint8')');
end