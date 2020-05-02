function [ rateMean, hh ] = globalRateTest(model, testSet, plotOption, clim)
%%%%predict rate from learned variational mean and variance of each user
% limitMax = 20;
% hh=[];
% for i=1:(limitMax +1)
%     for j=1:(i)
%         hh=[hh; [i-1,j-1]];
%     end;
% end;
hh = [];
for u=1:length(testSet);
    hh=[hh; testSet{u}];
end;
hh = unique(hh, 'rows');
if isfield(model, 'ratePar')
    ratePar = model.ratePar;
else
    ratePar.opt = 'pseudo';
    ratePar.cov = model.cov;
    ratePar.prior = model.prior;
    ratePar.ratio = model.ratio;
    ratePar.global.Mu = model.global.Mu;
    ratePar.GP.logtheta = model.GP.logtheta;
    ratePar.Xu = model.Xu;
    ratePar.maxRate = 3;
    ratePar.minRate = 0.01;
    ratePar.Kmm =  feval(ratePar.cov{:},ratePar.GP.logtheta, ratePar.Xu, []);
    ratePar.diagK = ratePar.Kmm(1);
    ratePar.KmmInv = ratePar.Kmm\eye(size(ratePar.Xu,1));
end;
if iscell(ratePar.cov{2}{1});
    ratePar.Knm = feval(ratePar.cov{2}{1}{:},ratePar.GP.logtheta(1:end-1),hh,ratePar.Xu);
else
    ratePar.Knm = feval(ratePar.cov{2}{1},ratePar.GP.logtheta(1:end-1),hh,ratePar.Xu);
end;
tmp = ratePar.prior.g + ratePar.Knm*ratePar.KmmInv*(ratePar.global.Mu - ratePar.prior.g);
rateMean = tmp.^2/ratePar.ratio;
if plotOption
    figure; surf_from_scatter(hh(:,1), hh(:,2), rateMean);
    if nargin>3
        caxis([clim(1), clim(2)]);
    end;
    xlabel('feature 1','FontSize',16);
    ylabel('feature 2','FontSize',16);
    set(gca,'FontSize',16);
end;

