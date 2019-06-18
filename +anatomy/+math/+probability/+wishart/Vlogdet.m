function out = Vlogdet(V, n, ~)
    K = size(V, 1);
    out = 0;
    for i=1:K
        out = out + psi(1, 0.5*(n+1-i));
    end
end