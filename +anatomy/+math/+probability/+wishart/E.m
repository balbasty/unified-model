function out = E(V, n, mode)
    if nargin < 3 || mode(1) ~= 'n'
        out = n*V;
    else
        out = V;
    end
end