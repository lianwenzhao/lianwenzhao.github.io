function [ x, y] = buildIntervalSet(eventsAll, M, Tstart, Tend, Tint, Tdur, Tstep )
%BUILDINTERVALSET Summary of this function goes here
U = length(eventsAll);
x = cell(1,U);
y = cell(1,U);
Lmax = max(Tint);
for u=1:U
    events = eventsAll{u};
    Nu = size(events,1);
    Tcur = Tstart+ Lmax;
    history = find((events(:,1)>= (Tcur-Lmax)).*(events(:,1) < Tcur)==1);
    count = max(history);
    while (Tcur+ Tdur < Tend)
        while ((~isempty(history)) & (events(history(1),1)< (Tcur-Lmax-1e-4)))
            history(1)=[];
        end;
        [~, feature] = featureExtractor(events(history, :), Tcur, Tint, M);
        x{u} = [x{u}; feature'];
        targetNum = 1;
        while (count+targetNum<= Nu)&(events(count+targetNum,1)<=Tcur+Tdur)
            targetNum = targetNum + 1;
        end;
        y{u} = [y{u}; targetNum-1];
        Tcur = Tcur + Tstep;
        if count == Nu
            break;
        end;
        while (events(count+1,1)<=Tcur)
%             disp(history)
%             disp(count)
            history = [history; count+1];
            count = count + 1;
            if count==Nu;
                break;
            end;
        end;
    end;
end
