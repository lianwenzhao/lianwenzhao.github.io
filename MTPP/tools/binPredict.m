function [inst, fea] = binPredict(model, events, Tint, Tstart, window)
M = 1;
Lmax= max(Tint);
label= zeros(M,1);
inst = [];
fea = [];
%%%%define rate
if isfield(model, 'ratePar')
    ratePar = model.ratePar;
    U = length(ratePar.Fm);
    fnProvided = 1;
else
    fnProvided = 0;
    U = length(model.events);
    ratePar.opt = 'pseudo';
    ratePar.cov = model.cov;
    ratePar.ratio = model.ratio;
    ratePar.prior = model.prior;
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
for u = 1:U
    event = events{u};
    idx = find(event(:,1)>=Tstart, 1, 'first');
    Tcur = event(idx,1);
    while  (Tcur<max(event(:,1))-window)
        idx= find((event(:,1)>Tcur).*(event(:,1)<=Tcur+window)==1);
        for m=1:M
            label(m) = sum(event(idx,2)==m)>0;
        end;
        hist= find((event(:,1)<Tcur).*(event(:,1)>=Tcur-Lmax)==1);
        [duration, feature] = featureExtractor(event(hist, :), Tcur, Tint, M);
        [rate] = feature2rate(feature, ratePar, u, fnProvided);
        fea = [fea; feature'];
        cumHazard = 0;
        Tnext= Tcur + window;
        while (Tcur+ duration < Tnext)
            cumHazard = cumHazard + duration* rate;
            Tcur = Tcur + duration;
            [duration, feature] = featureExtractor(event(hist, :), Tcur, Tint, M);
            [rate] = feature2rate(feature, ratePar, u, fnProvided);
        end
        cumHazard = cumHazard + rate * (Tnext - Tcur);
        Tcur = Tnext;
        inst = [inst; label, 1-exp(-cumHazard)];
    end;
end;
