function type = range2int(maxval,minval)
%__________________________________________________________________________
% anatomy.segment.gmm.lib.range2int
%--------------------------------------------------------------------------
% FORMAT type = gmm.lib.range2int(maxval, minval)
%
% Find the best suited integer type to store integer values in the range
% [minval,maxval]
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging


if nargin < 2
    minval = 0;
end

type     = 'int';
unsigned = minval >= 0;
if unsigned
    type   = ['u' type];
    minval = 0;
else
    minval = numel(dec2base(-minval,2));
end
maxval = numel(dec2base(maxval,2));
nbits  = max(minval,maxval);
if unsigned
    nbits = nbits + 1;
end
if nbits <= 8
    type = [type '8'];
elseif nbits <= 16
    type = [type '16'];
elseif nbits <= 32
    type = [type '32'];
elseif nbits <= 64
    type = [type '64'];
else
    type = 'double';
end