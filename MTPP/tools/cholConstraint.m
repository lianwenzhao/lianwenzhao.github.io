function [AA, bb, AAe, bbe] = cholConstraint(M)
tmp = zeros(M,M);
if nargout==2
    cnt = 1;
    for i=1:M
        for j=i:M
            tmp(j, i) = cnt;
            cnt = cnt+1;
        end;
    end;
    AA = zeros(M, M+M*(M+1)/2);
    for i=1:M
        AA(i, M + tmp(i,i)) = -1;
    end;
    bb = -1e-2*ones(M, 1);
else
    cnt = 1;
    for i=1:M
        for j=1:M
            tmp(j, i) = cnt;
            cnt = cnt+1;
        end;
    end;
    cnt1 = 1;
    cnt2 = 1;
    AA = zeros(M, M + M^2);
    AAe = zeros(M*(M-1)/2, M + M^2);
    for i=1:M
        for j = 1:M
            if (i==j)
                AA(cnt1, M+tmp(j, i)) = -1;
                cnt1 = cnt1+1;
            elseif (i>j)
                AAe(cnt2, M+tmp(j,i)) = 1;
                cnt2 = cnt2+1;
            end;
        end;
    end;
    bb = -1e-2*ones(M, 1);
    bbe = zeros(M*(M-1)/2, 1);
    
end;