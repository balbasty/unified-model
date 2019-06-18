function C = mtimes(A, B, dims)
%__________________________________________________________________________
% anatomy.math.fields.mtimes
%--------------------------------------------------------------------------
% FORMAT C = mtimes(A, B, [dims])
% A    - A field of matrices of dimensions [... M N]
% B    - A field of matrices of dimensions [... N P]
% C    - A field of matrices of dimensions [... M P]
% dims - Dimensions of the nd-arrays corresponding to the matrix dimensions
%           default: [2 3]
%
% Matrix multiplication between fields of matrices, with implicit
% expansion enabled.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

    % By default: assume a three-dimensional field
    if nargin < 3 || isempty(dims)
        dims = [2 3];
    end
    % If only one dimension provided, assum we permute with the next one
    if numel(dims) == 1
        dims = [dims dims+1];
    end
    % If more than two dimensions provided, issue a warning
    if numel(dims) > 2
        warning('Only two dimensions required. I will discard the remaining ones')
        dims = dims(1:2);
    end
    
    % Check that dimensions match correctly
    dima = [size(A) 1 1];
    dimb = [size(B) 1 1];
    if dima(dims(2)) ~= dimb(dims(1))
        error('A and B do not have comformable sizes (%d ~= %d)', ...
            dima(dims(2)), dimb(dims(1)));
    end
    dimmata = dima(dims);
    dimmatb = dimb(dims);
    dima(dims) = 1;
    dimb(dims) = 1;
    ok = all((dima == dimb) | (dima == 1) | (dimb == 1));
    if ~ok
        error('The corresponding (field) dimensions of A and B must be equal to each other, or equal to one')
    end
    
    dim  = max(dima, dimb);
    dimc = dim;
    dimc(dims) = [dimmata(1) dimmatb(2)];
    C = zeros(dimc, 'like', A);
    
    for i=1:dimmata(1)
        for j=1:dimmatb(2)
            % Initialise C(:,i,j)
            Cij = zeros(1, 'like', C);
            for k=dimmata(2)
                % Read A(:,i,k)
                Sa = struct;
                Sa.type = '()';
                Sa.subs = repmat({':'}, [1 numel(dima)]);
                Sa.subs(dims) = [i k];
                Aik = subsref(A, Sa);
                % Read B(:,k,j)
                Sb = struct;
                Sb.type = '()';
                Sb.subs = repmat({':'}, [1 numel(dimb)]);
                Sb.subs(dims) = [k j];
                Bkj = subsref(B, Sb);
                % Compute product and accumulate
                Cij = Cij + bsxfun(@times, Aik, Bkj);
                clear Aik Bkj
            end
            % Store result in C
            Sc = struct;
            Sc.type = '()';
            Sc.subs = repmat({':'}, [1 numel(dimc)]);
            Sc.subs(dims) = [i j];
            C = subsasgn(C, Sc, Cij);
            clear Cij
        end
    end
    
end