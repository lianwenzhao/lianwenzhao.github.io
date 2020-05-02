function [ pieceAll ] = stream2instance(eventAll, Tint, Tend, M)
%STREAM2INSTANCE Summary of this function goes here
if iscell(eventAll)
    U = length(eventAll);
    pieceAll = cell(1, U);
    for u=1:U
        if isempty(eventAll{u})
            pieces.feature = [];
            pieces.label = [];
            pieces.duration = [];
            pieces.time = [];
        else
            events = eventAll{u};
            Lmax =  max(Tint);
            Tcur = Lmax;
            while ((Tcur+Lmax)<events(1,1))
                Tcur = Tcur + Lmax;
            end;
            first= find(events(:,1)>Tcur, 1, 'first');%%first event starts after Lmax instead of >=
            lastVisit = first-1;
            hist = [1: lastVisit];
            instances=[];
            stamp=[];
            while (Tcur<events(end,1))
                %     if (lastVisit==size(events,1))
                %         disp(num2str(Tcur));
                %     end;
                [duration, feature] = featureExtractor(events(hist,:), Tcur, Tint, M);
                if Tcur+duration>events(lastVisit+1,1)
                    lastVisit = lastVisit+1;
                    instances=[instances; events(lastVisit,2), events(lastVisit,1)-Tcur, feature' ];
                    hist= [hist lastVisit];
                    Tcur= events(lastVisit,1);
                    stamp=[stamp; Tcur];
                    if (lastVisit < size(events,1))
                        while (events(lastVisit+1,1)==Tcur)
                            lastVisit = lastVisit+1;
                            instances=[instances; events(lastVisit,2), events(lastVisit,1)-Tcur, feature' ];
                            stamp=[stamp; Tcur];
                            hist= [hist lastVisit];
                        end;
                    end;
                else
                    Tcur= Tcur+ duration;
                    instances=[instances; 0, duration, feature' ];
                    stamp=[stamp; Tcur];
                end;
                while ((~isempty(hist)) && (events(hist(1),1)< (Tcur-Lmax-1e-4)))
                    hist(1)=[];
                end;
            end;
            while (Tcur<Tend)
                [duration, feature] = featureExtractor(events(hist,:), Tcur, Tint, M);
                Tcur= Tcur+ duration;
                instances=[instances; 0, duration, feature' ];
                stamp=[stamp; Tcur];
                while ((~isempty(hist)) && (events(hist(1),1)< (Tcur-Lmax-1e-4)))
                    hist(1)=[];
                end;
            end;
            if (stamp(end)>Tend)
                stamp(end)= Tend;
                instances(end,:)= [0, Tend-stamp(end-1), instances(end,3:end)];
            end;
            pieces.feature = instances(:, 3:end);
            pieces.label = instances(:,1);
            pieces.duration = instances(:,2);
            pieces.time = stamp;
        end;
        pieceAll{u}= pieces;
    end;
else
    if isempty(eventAll)
        pieces.feature = [];
        pieces.label = [];
        pieces.duration = [];
        pieces.time = [];
    else
        events = eventAll;
        Lmax =  max(Tint);
        Tcur = Lmax;
        while ((Tcur+Lmax)<events(1,1))
            Tcur = Tcur + Lmax;
        end;
        first= find(events(:,1)>=Tcur, 1, 'first');
        lastVisit = first-1;
        hist = [1: lastVisit];
        instances=[];
        stamp=[];
        while (Tcur<events(end,1))
            %     if (lastVisit==size(events,1))
            %         disp(num2str(Tcur));
            %     end;
            [duration, feature] = featureExtractor(events(hist,:), Tcur, Tint, M);
            if Tcur+duration>events(lastVisit+1,1)
                lastVisit = lastVisit+1;
                instances=[instances; events(lastVisit,2), events(lastVisit,1)-Tcur, feature' ];
                hist= [hist lastVisit];
                Tcur= events(lastVisit,1);
                stamp=[stamp; Tcur];
                if (lastVisit < size(events,1))
                    while (events(lastVisit+1,1)==Tcur)
                        lastVisit = lastVisit+1;
                        instances=[instances; events(lastVisit,2), events(lastVisit,1)-Tcur, feature' ];
                        stamp=[stamp; Tcur];
                        hist= [hist lastVisit];
                    end;
                end;
            else
                Tcur= Tcur+ duration;
                instances=[instances; 0, duration, feature' ];
                stamp=[stamp; Tcur];
            end;
            while ((~isempty(hist)) && (events(hist(1),1)< (Tcur-Lmax-1e-4)))
                hist(1)=[];
            end;
        end;
        while (Tcur<Tend)
            [duration, feature] = featureExtractor(events(hist,:), Tcur, Tint, M);
            Tcur= Tcur+ duration;
            instances=[instances; 0, duration, feature' ];
            stamp=[stamp; Tcur];
            while ((~isempty(hist)) && (events(hist(1),1)< (Tcur-Lmax-1e-4)))
                hist(1)=[];
            end;
        end;
        if (stamp(end)>Tend)
            stamp(end)= Tend;
            instances(end,:)= [0, Tend-stamp(end-1), instances(end,3:end)];
        end;
        pieces.feature = instances(:, 3:end);
        pieces.label = instances(:,1);
        pieces.duration = instances(:,2);
        pieces.time = stamp;
    end;
    pieceAll = pieces;
end;
