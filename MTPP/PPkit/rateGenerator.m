function [ ratePar ] = rateGenerator(rateOption, U)
if nargin<2
    U = 1;
end;
%RATEGENERATOR Summary of this function goes here
if strcmp(rateOption,'linear')
    ratePar.opt = 'linear';
    ratePar.weights = [0.25 -0.1]';
    ratePar.weights = [ratePar.weights [0.1; 0.3]];
    ratePar.scaling = 2;
    ratePar.maxRate = 1;
    ratePar.minRate = 0.1;
elseif strcmp(rateOption, 'sigmoid')
    ratePar.opt = 'sigmoid';
    ratePar.weights = [0.25 -0.1]';
    ratePar.scaling = 1;
    ratePar.maxRate = 1;
    ratePar.minRate = 0.1;
elseif strcmp(rateOption, 'pseudo')
    ratePar.opt = 'pseudo';
    %     ratePar.cov = {'covSum',{'covSEard','covNoise'}};
    covM3 = {@covMaternard,3};
    ratePar.cov = {'covSum',{covM3,'covNoise'}};
    ss = 0;
    mm = 2;
    ll = 9;
    bb = 7;
    SS = 0;
    MM = 4;
    LL = 10;
    BB = 15;
    ratePar.Xu = [SS ss; MM ss; MM mm;
        LL ss; LL mm; LL ll];
    ratePar.prior.g = 3;
    ratePar.ratio = 10;
    ratePar.GP.logtheta = log([4, 4, 1.2, 0.01]);
    ratePar.Kmm =  feval(ratePar.cov{:},ratePar.GP.logtheta, ratePar.Xu, []);
    ratePar.diagK = ratePar.Kmm(1);
    ratePar.KmmInv = ratePar.Kmm\eye(size(ratePar.Xu,1));

    ratePar.global.Mu = sqrt([0.1 1.22 1.22...
            1.22 1.22 0.1]' * ratePar.ratio);
    for u=1:U
        ratePar.Fm{u} = sqrt([0.1 1.22 1.22...
            1.22 1.22 0.1]' * ratePar.ratio);
        ratePar.KinvF{u} = ratePar.KmmInv * (ratePar.Fm{u}-ratePar.prior.g);
    end;
    ratePar.maxRate = 3;
    ratePar.minRate = 0.1;
    limitMax = 23;
    hh=[];
    for i=1:(limitMax +1)
        for j=1:(i)
            hh=[hh; [i-1,j-1]];
        end;
    end;
    ratePar.Knn =  feval(ratePar.cov{:},ratePar.GP.logtheta, hh, []);
    if iscell(ratePar.cov{2}{1});
        ratePar.Knm = feval(ratePar.cov{2}{1}{:},ratePar.GP.logtheta(1:end-1),hh,ratePar.Xu);
    else
        ratePar.Knm = feval(ratePar.cov{2}{1},ratePar.GP.logtheta(1:end-1),hh,ratePar.Xu);
    end;
    tmpVar = ratePar.Knn - ratePar.Knm*ratePar.KmmInv*ratePar.Knm';
    cholVar = chol(tmpVar)';
    for u=1:U
        fn = 1e3;
        while ((max(fn.^2/ratePar.ratio)>ratePar.maxRate)||min(fn)<0)
            fn = (ratePar.prior.g + ratePar.Knm*ratePar.KinvF{u} + cholVar*randn(size(hh,1),1));
        end;
        ratePar.fnTable{u} = containers.Map;
        cnt = 0;
        for i=1:(limitMax +1)
            for j=1:(i)
                cnt= cnt+1;
                ratePar.fnTable{u}(num2str(hh(cnt,:))) = fn(cnt);
            end;
        end;
    end;
end

