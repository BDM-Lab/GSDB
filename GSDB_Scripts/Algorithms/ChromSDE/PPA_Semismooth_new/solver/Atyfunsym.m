%%*********************************************************
%% Atyfun: compute sum_{k=1}^m yk*Ak.
%%
%%  Q = Atyfun(blk,At,y);
%%**********************************************************

function Q = Atyfunsym(blk,At,y,isspAy)

if (nargin < 4); isspAy = ones(size(blk,1),1); end
Q = cell(size(blk,1),1);
%%
if iscell(At)
    for p = 1:size(blk,1)
        pblk = blk(p,:);
        if strcmp(pblk{1},'s')
            Q{p} = mexsmat(pblk,At{p,1}*y);
        elseif strcmp(pblk{1},'l')
            Q{p} = At{p,1}*y;
        end
    end
else
    if strcmp(blk{1,1},'s')
        Q = mexsmat(blk,At*y);
    elseif strcmp(blk{1,1},'l')
        Q = At*y;
    end
end
%%*********************************************************

