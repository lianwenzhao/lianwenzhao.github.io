function llPred = predLikelihood(rates, stamps, events, Tstart, Tend)
%%%comput likelihood for piecewise constant point process
U = length(events);
llPred = 0;
for u=1:U
    idx = find(stamps{u}>=Tstart,1,'first');
    stamps{u} = stamps{u}(idx:end,:);
    rates{u} = rates{u}(idx:end,:);
    N = size(events{u},1);
    C = length(stamps{u});
    j = 1;
    i = 1;
    while (events{u}(i)<stamps{u}(2))
        llPred = llPred - rates{u}(1);
        i = i+1;
    end;
    flag = 1;
    while (flag)
        while ((stamps{u}(j+1)-events{u}(i,1))<1e-6) 
            llPred = llPred - (stamps{u}(j+1) - stamps{u}(j)) * rates{u}(j);
            j = j+1;
            if j>=C
                flag= 0;
                break;
            end;
        end;
        if (abs(stamps{u}(j)-events{u}(i,1))>1e-6) %(stamps{u}(j)~= events{u}(i))
            disp('error!');
        else
            llPred = llPred + log(rates{u}(j-1)+eps);
            i = i+1;
            if j==C
                break;
            end;
            if i>N
                while (stamps{u}(j+1) <=stamps{u}(C))
                    llPred = llPred - (stamps{u}(j+1) - stamps{u}(j)) * rates{u}(j);
                    j = j+1;
                    if j==C
                        flag= 0;
                        break;
                    end;
                end;
            end;
        end;
    end;
    
end

