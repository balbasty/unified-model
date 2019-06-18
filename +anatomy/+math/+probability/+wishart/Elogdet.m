function out = Elogdet(V, n, mode)
    
    import anatomy.math.matrix                 % posdef.logdet
    import anatomy.math.probability            % digamma
    import anatomy.math.probability.wishart    % Elogdet

    if nargin < 3 || mode(1) ~= 'n'
        K   = size(V, 1);
        out = digamma(0.5*n, K) + K*log(2) + posdef.logdet(V);
    else
        out = Elogdet(V/n, n);
    end
end