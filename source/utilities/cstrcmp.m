function cmp = cstrcmp(a, b)
%CSTRCMP  Compare two character arrays according to ASCII dictionary order in C
%style.
%
%   INPUTS:
%
%   -|a|, |b|: character arrays, not necessarily of the same length
%
%   OUTPUT:
%
%   -|cmp|: integer holding the result of the comparison. Equal to 0 when the
%    two inputs are equal, positive when |a| is greater than |b|, and negative
%    when |a| is less than |b|.

% Force the arrays to equal length.
x = char({a;b});

% Subtract one from the other.
d = x(1,:) - x(2,:);


if any(d)
    cmp = d(find(d, 1));
else
    cmp = 0;
end
    
end