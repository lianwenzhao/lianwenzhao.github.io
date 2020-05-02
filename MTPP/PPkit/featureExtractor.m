function [duration, feature] = featureExtractor(events, Tcur, Tint, M)
R= length(Tint);
duration = max(Tint);
if (isempty(events))
    feature= zeros(M*R,1);
else
    lb= Tcur- Tint;%%lower bound of the feature intervals
    featurePerMark=cell(M,1);
    for m=1:M
        featurePerMark{m}=zeros(R,1);
    end;
    for i=1:size(events, 1)
        tmp = (events(i,1)-lb);
        featurePerMark{events(i,2)} = featurePerMark{events(i,2)}+ (tmp>1e-5);
        duration = min([tmp(tmp>1e-5); duration]);
    end;
%     if (duration == 1e10)
%         disp('error!');
%         %     pause;
%     end;
    feature = cell2mat(featurePerMark);
end;


