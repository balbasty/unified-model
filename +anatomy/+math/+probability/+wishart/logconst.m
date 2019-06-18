function out = logconst(V,n)
    
    import anatomy.math.matrix

    M     = size(V,1);
    if M >= n+1
       %warning('SPM:Wishart','Can not normalise a Wishart distribution (M=%d, nu=%f)', M,nu);
        out = 0;
        return;
    end
    lGamM = M*(M-1)/4*log(pi);
    for m=1:M
        lGamM = lGamM + gammaln((n + 1 - m)/2);
    end
    out = 0.5*n*(posdef.logdet(V) + M*log(2)) + lGamM;
end