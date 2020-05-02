function [ rateMean] = varPpPredict(model, testing)
%%%%predict rate from learned variational mean and variance of each user
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
rateMean = cell(1,U);
for u=1:U
    for i=1:size(testing{u},1)
        [tmp] = feature2rate(testing{u}(i,:)', ratePar, u, fnProvided);
        rateMean{u} = [rateMean{u}; tmp];
    end;
end;

