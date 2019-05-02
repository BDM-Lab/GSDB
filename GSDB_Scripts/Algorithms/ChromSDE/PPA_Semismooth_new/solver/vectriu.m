%%***********************************************************
%% vectriu: stack the upper triangular part of a matrix
%%          into a vector.
%%
%% xvec = vectriu(blk,x);
%%
%%***********************************************************

function xvec = vectriu(blk,x)

if ~iscell(x)
    r = blk{1,2};
    ismt = isempty(x);
    if strcmp(blk{1,1},'s') & (r > 0);
        r2 = r*(r+1)/2;
        xvec = zeros(r2,1);
        idx = 0;
        for j = 1:r
            xvec(idx+[1:j]) = x(1:j,j);
            idx = idx + j;
        end
    elseif strcmp(blk{p,1},'l') & (r > 0);
        xvec = x;
    end
else
    xvec = [];
    for p = 1:size(blk,1);
        xvec = [xvec; vectriu(blk(p,:),x{p})];
    end
end
%%***********************************************************
