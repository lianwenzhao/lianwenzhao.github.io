function [ rateMean, count] = ratePredict(model, testing)
%%%%predict rate from learned variational mean and variance of each user
if isfield(model, 'ratePar')
    U = length(model.ratePar.fnTable);
    rateMean=cell(1,U);
    rateVar = cell(1,U);
    for u = 1:U
        Xnu = unique(testing{u}, 'rows');
        for i=1:size(Xnu,1)
            rateMean{u} = [rateMean{u}; (model.ratePar.fnTable{u}(num2str(Xnu(i,:)))).^2/model.ratePar.ratio];
        end;
    end;
else
    pseudo.x = model.Xu;
    pseudo.y = model.var.Mu;
    model.Kmm = feval(model.cov{:}, model.GP.logtheta, pseudo.x);
    U = length(model.Nu);
    if nargin<2
        for u = 1:U
            Xnu = unique(model.events{u}.feature, 'rows');
            model.Knm{u} = feval(mid, model.GP.logtheta(1:end-1), Xnu, pseudo.x);
            model.Knn{u} = feval(model.cov{:}, model.GP.logtheta, Xnu);
        end;
    else
        for u = 1:U
            Xnu = unique(testing{u}, 'rows');
            if iscell(model.cov{2}{1});
                model.Knm{u} = feval(model.cov{2}{1}{:}, model.GP.logtheta(1:end-1), Xnu, pseudo.x);
            else
                model.Knm{u} = feval(model.cov{2}{1}, model.GP.logtheta(1:end-1), Xnu, pseudo.x);
            end;
            model.Knn{u} = feval(model.cov{:}, model.GP.logtheta, Xnu);
        end;
    end;
    rateMean=cell(1,U);
    for u = 1:U
        tmp = model.Knm{u}/ model.Kmm;
        fMean = model.prior.g + tmp * (pseudo.y{u}-model.prior.g);
        fVar = diag(model.Knn{u} - tmp* model.Knm{u}');
        rateMean{u} = (fMean).^2/model.ratio;%rateMean{u} = exp(fMean+ 0.5*fVar);
    end;
end;
%%%count number of appearance of unique features 
count = cell(1,U);
for u = 1:U
    Xnu = unique(testing{u}, 'rows');
    tmpTable = containers.Map;
    for i=1:size(Xnu,1)
        tmpTable(num2str(Xnu(i,:))) = 0;
    end;
    for j=1:size(testing{u},1)
        tmpTable(num2str(testing{u}(j,:))) = tmpTable(num2str(testing{u}(j,:))) + 1;
    end;
    count{u} = cell2mat(tmpTable.values');
end;
