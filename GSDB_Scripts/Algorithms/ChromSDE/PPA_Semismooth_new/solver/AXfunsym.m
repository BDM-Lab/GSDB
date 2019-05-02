%%*********************************************************
%% AXfun: compute AX(k) = <Ak,X>, k = 1:m
%%
%% AX = AXfun(blk,At,X);
%%
%%**********************************************************

function AX = AXfunsym(blk,At,X)

if iscell(At)
    m = size(At{1},2);
    AX = zeros(m,1);
    for p = 1:size(blk,1);
        pblk = blk(p,:);
        if strcmp(pblk{1},'s')
            if (length(pblk{2}) == 1)
                AX = AX + (mexsvec(pblk,X{p})'*At{p,1})';
            else
                AX = AX + (mexsvec(pblk,sparse(X{p}))'*At{p,1})';
            end
        elseif strcmp(pblk{1},'l')
            AX = AX + (X{p}'*At{p,1})';
        end
    end
else
    if strcmp(blk{1,1},'s')
        if (length(blk{1,2})==1)
            AX = (mexsvec(blk,X)'*At)';
        else
            AX = (mexsvec(blk,sparse(X))'*At)';
        end
    elseif strcmp(blk{1,1},'l')
        AX = (X'*At)';
    end
end
%%*********************************************************
