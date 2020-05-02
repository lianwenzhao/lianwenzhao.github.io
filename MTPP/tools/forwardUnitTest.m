function [sse, probHit, features, trueCount] = forwardUnitTest(Tstart, Tend, Tint, M, model, eventsObs, MCiter, plotOption)
%%%forward sampling given the model and history
if isfield(model, 'ratePar')
    ratePar = model.ratePar;
    U = length(ratePar.Fm);
    fnProvided = 1;
else
    fnProvided = 0;
    U = length(model.events);
    ratePar.opt = 'pseudo';
    ratePar.cov = model.cov;
    ratePar.prior = model.prior;
    ratePar.ratio = model.ratio;
    ratePar.GP.logtheta = model.GP.logtheta;
    ratePar.Xu = model.Xu;
    ratePar.maxRate = 3;
    ratePar.minRate = 0.01;
    ratePar.Kmm =  feval(ratePar.cov{:},ratePar.GP.logtheta, ratePar.Xu, []);
    ratePar.diagK = ratePar.Kmm(1);
    ratePar.KmmInv = ratePar.Kmm\eye(size(ratePar.Xu,1));
    for u=1:U
        ratePar.Fm{u} = model.var.Mu{u}-model.prior.g;% + 0.5*diag(model.var.Sigma{u}));
        ratePar.KinvF{u} = ratePar.KmmInv * ratePar.Fm{u};
        ratePar.SigmaInv{u} = model.var.Sigma{u}\eye(model.m);
    end;
end;
eventAll = cell(MCiter, U);
instAll = cell(MCiter,U);
rateAll = cell(MCiter, U);
eventCount = zeros(MCiter, U);
if strcmp(plotOption, 'on');    fig = figure;end;
Lmax= max(Tint);
eventsInit = cell(1,U);
trueCount = zeros(1,U);
features = zeros(length(Tint), U);
for user = 1:U
    eventsInit{user} = eventsObs{user}(((eventsObs{user}(:,1)>= (Tstart-Lmax)).*(eventsObs{user}(:,1) < Tstart)==1),:);
    trueCount(user) = sum((eventsObs{user}(:,1)>=Tstart).*(eventsObs{user}(:,1)<Tend));
    [~, features(:,user)] = featureExtractor(eventsInit{user}, Tstart, Tint, M);
end;
for mc = 1:MCiter
    for user=1:U
        Tcur = Tstart;
        events= []; %%time stamp; mark type
        m=1;
        events = eventsInit{user};%%dynamically store the index of events for history
        count0 = size(events, 1);
        count = count0;
        history = [1:1:count];
        inst = [];
        ratePiece = [];
        while  (Tcur<Tend)
            u = rand(M,1);
            while (Tcur<Tend)
                while ((~isempty(history)) && (events(history(1),1)< (Tcur-Lmax-1e-4)))
                    history(1)=[];
                end;
                [duration, feature] = featureExtractor(events(history, :), Tcur, Tint, M);
                [rate] = feature2rate(feature, ratePar, user, fnProvided);
                if ~isreal(rate)
                    disp('wrong');
                end;
                %             [rate] = feature2rate(feature, ratePar, clu(user));
                uPiece= 1- exp(- rate* duration);
                if ((sum(uPiece > u)>0))%%generate an event
                    break;
                else%%%no event
                    Tcur = Tcur+ duration;
                    u = (u-uPiece)./(1- uPiece);
                    inst= [inst; 0, duration, feature'];
                    ratePiece = [ratePiece; rate'];
                end;
            end;
            [~, m] = max(uPiece - u);%find the mark type
            Tcur = Tcur - log( eps+1- u(m))/(eps+rate(m));
            inst= [inst; m,  -log( eps+1- u(m))/(eps+rate(m)), feature'];
            ratePiece = [ratePiece; rate'];
            if Tcur > Tend
                break;
            end;
            events = [events; [Tcur, m]];
            count = count+1;
            history = [history count];
        end
        assert(~isempty(inst), 'error in generating events!');
        idx = find((Tstart+cumsum(inst(:,2))<Tend), 1, 'last');
        if (length(idx)==0)
            idx = 0;
        end;
        tmp1 = inst(idx+1,:);
        tmp2 = ratePiece(idx+1,:);
        inst= inst(1:idx,:);
        ratePiece = ratePiece(1:idx, :);
        inst = [inst; 0, Tend - sum(inst(:,2))- Tstart, tmp1(3:end)];
        ratePiece = [ratePiece; tmp2];
        ratePiece = [Tstart+ [0; cumsum(inst(:,2))] [ratePiece(1,:); ratePiece]];
        eventAll{mc, user}=events(count0+1:end,:);
        instAll{mc, user} = inst;
        rateAll{mc, user} = ratePiece;
        eventCount(mc, user) = size(eventAll{mc, user}, 1);
    end;
end;
probHit = zeros(1,U);
sse = zeros(1, U);
for user = 1:U
    cc = unique(eventCount(:,user));
    nn= zeros(1,length(cc));
    for j=1:length(cc);
        nn(j) = sum(cc(j)==eventCount(:,user));
    end;
%     [nn, cc] = hist([eventCount(:,user)], 10);
    pp = nn/MCiter;
    sse(user) = sum(pp*((cc-trueCount(user)).^2));
    [~, idx ]= min(abs(cc-trueCount(user)));
    probHit(user) = pp(idx);
end;
for user = 1:3
    if strcmp(plotOption,'on')
        for m=1:M
            subplot(M*3,1,m+M*(user-1));
            tmp=eventAll{1,user}(eventAll{1,user}(:,2)==m,1);
            stem(tmp, ones(length(tmp),1),'o');
            axis([Tstart, Tend, 0, 1.5]);
        end;
        plotPiecewiseRate(rateAll{mc, user} , fig, 'r');
    end;
end;
