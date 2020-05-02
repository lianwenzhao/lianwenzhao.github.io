function [ rmse, probHit, features, trueCount ] = forwardTest( TtestStart, TtestDur, TtestEnd, Tint, M, model, eventsAll, maxNumTest, MCiter)
%FORWARDTEST Summary of this function goes here
if ~exist('MCiter','var'),  MCiter = 100; end
U = length(eventsAll);
sse = zeros(1,U);
probHit =[];
features = cell(1,U);
trueCount = [];
TpredStart = TtestStart;
count = 0;
for i =1:maxNumTest
    TpredEnd = TpredStart + TtestDur; %TtestEnd;
    if TpredEnd > TtestEnd;break;end;
    count = count + 1;
    [sseSingle, probHitSingle, featuresSingle, trueCountSingle] = ...
    forwardUnitTest(TpredStart, TpredEnd, Tint, M, model, eventsAll, MCiter, 'off');
    for u=1:U
        features{u} = [features{u}; featuresSingle(:,u)'];
    end;
    trueCount = [trueCount; trueCountSingle];
    sse = sse + sseSingle;
    probHit = [probHit; probHitSingle];
    TpredStart = TpredStart + TtestDur;
    disp(['Testing: Tcur ' num2str(TpredStart) ' Tend ' num2str(TtestEnd)])
end
rmse = sqrt(sse/count);

