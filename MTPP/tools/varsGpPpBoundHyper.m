function  [f, df] = varsGpPpBoundHyper(W, model, flag)
%%%%compute f and df of lower bound w.r.t. theta
%%%%theta includes inducing variables, GP hyper-parameters

% PRELIMINARIES
if nargin == 2
    flag = 1;
end

% place W to the model structure
model = returnOptimizedParams(model, W, flag);

U = length(model.events);

% COVARIANCE FUNCTION QUANTITIES needed in the sparse method: K_mm, Knm, diag(Knn)
%
% inducing variable parameters are pseudo-inputs
if strcmp(model.indType, 'pseudoIns')
    model.Kmm = feval(model.cov{:}, model.GP.logtheta,model.Xu);
    if iscell(model.cov{2}{1});
        for u=1:U
            model.Knm{u} = feval(model.cov{2}{1}{:},model.GP.logtheta(1:end-1),model.events{u}.feature,model.Xu);
        end;
    else
        for u=1:U
            model.Knm{u} = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.events{u}.feature,model.Xu);
        end;
    end;
end

model.diagKnn = sum(exp(2*model.GP.logtheta(end-1:end)) );%%%simply the sigma_f^2
if model.GP.constDiag == 1
    for u=1:U
        model.TrKnn{u} = size(model.events{u}.feature,2)*model.diagKnn;
    end;
else
    error('Knn not properly defined! ');
end

for u = 1:U
    w{u} = model.events{u}.label;
    v{u} = model.events{u}.duration;
end;
invKmm = model.Kmm\eye(model.m);
f = 0;
for u = 1:U
    KnmInvKmm{u} = model.Knm{u}*invKmm;
    meanRate{u} = KnmInvKmm{u}*(model.var.Mu{u}-model.prior.g) + model.prior.g;
    meanRate2{u} = meanRate{u}.^2;
    varRate{u} = (1+1/model.prior.xi)*(model.diagKnn - diag(KnmInvKmm{u}*model.Knm{u}'))...
        + diag(KnmInvKmm{u}*model.var.Sigma{u}*KnmInvKmm{u}');
    assert(min(varRate{u})>0, 'not positive definate!');
    for i=1:model.Nu(u)
        if (-meanRate2{u}(i)/2/varRate{u}(i)>0)
            disp('wrong!');
        end;
        [gg, ~] = queryGz(-meanRate2{u}(i)/2/varRate{u}(i), model.gz, model.dgz);
        if w{u}(i)
            f = f -  gg+ log(0.5*varRate{u}(i));
        end;
    end;
    f = f - v{u}'*(meanRate2{u}+varRate{u});
    f = f - 0.5*model.prior.xi* ((model.global.Mu-model.prior.g)'*invKmm*(model.global.Mu-model.prior.g))...
        - 0.5*trace(invKmm* (model.var.Sigma{u} + (model.var.Mu{u}-model.global.Mu) * (model.var.Mu{u}-model.global.Mu)' ) );
end;
f = f  - (U+1)/2 * logdet(model.Kmm); 
f = -f;

df = zeros(length(W),1);%zeros(model.D*model.m+model.GP.nParams+1,1);
if (flag==3)%update GP hyper parameters only
    numPar = model.GP.nParams;
    for k=1:(model.GP.nParams)
        DkmmDx{k} = feval(model.cov{:}, model.GP.logtheta, model.Xu, [], k);%%noise term is not considered@covSum, @covNoise
        for u=1:U
            if (k<model.GP.nParams)
                if iscell(model.cov{2}{1})
                    DknmDx{k}{u} = feval(model.cov{2}{1}{:}, model.GP.logtheta(1:end-1), model.events{u}.feature,model.Xu,  k);
                else
                    DknmDx{k}{u} = feval(model.cov{2}{1}, model.GP.logtheta(1:end-1), model.events{u}.feature,model.Xu,  k);                    
                end;
            else
                DknmDx{k}{u} = zeros(model.Nu(u), model.m);
            end;
            DknnDx{k}{u} = feval(model.cov{:}, model.GP.logtheta, model.events{u}.feature, [],  k);
        end;
    end;
    
    % update both inducing variables and GP hyper parameters
else
    if (flag==2)
        numPar = model.D*model.m;
    elseif (flag ==1)
        numPar = (model.D*model.m+model.GP.nParams);
    end;
    for i = 1:model.m
        tmp= - repmat(model.Kmm(i,:), [model.D,1]).*...%%only the i-th row and i-th col has non-zeros
            (diag(exp(-model.GP.logtheta(1:model.D)*2))*(repmat(model.Xu(i,:)',[1,model.m]) - model.Xu'));
        for d=1:model.D
            DkmmDx{(i-1)*model.D+d} = zeros(model.m, model.m);
            DkmmDx{(i-1)*model.D+d}(:,i) = tmp(d,:)';
            DkmmDx{(i-1)*model.D+d}(i,:) = tmp(d,:);
        end;
        for u=1:U
            tmp = - repmat(model.Knm{u}(:,i), [1, model.D]).*...%%only the i-th col has non-zeros
                ((repmat(model.Xu(i,:),[model.Nu(u), 1]) - model.events{u}.feature)*diag(exp(-model.GP.logtheta(1:model.D)*2)));
            for d=1:model.D
                DknmDx{(i-1)*model.D+d}{u} = zeros(model.Nu(u), model.m);
                DknmDx{(i-1)*model.D+d}{u}(:,i) = tmp(:,d);
                DknnDx{(i-1)*model.D+d}{u} = 0;
            end;
        end;
    end;
    for k=(model.D*model.m+1):numPar
        DkmmDx{k} = feval(model.cov{:}, model.GP.logtheta, model.Xu, [], k-model.D*model.m);
        for u=1:U
            if (k<model.D*model.m+model.GP.nParams)
                if iscell(model.cov{2}{1})
                    DknmDx{k}{u} = feval(model.cov{2}{1}{:}, model.GP.logtheta(1:end-1), model.events{u}.feature,model.Xu, k-model.D*model.m);
                else
                    DknmDx{k}{u} = feval(model.cov{2}{1}, model.GP.logtheta(1:end-1), model.events{u}.feature,model.Xu, k-model.D*model.m);
                end;
            else
                DknmDx{k}{u} = zeros(model.Nu(u), model.m);
            end;
            DknnDx{k}{u} = feval(model.cov{:}, model.GP.logtheta, model.events{u}.feature, [],  k-model.D*model.m);
        end;
    end;
end;

for k = 1:numPar
    DkmmInvDx = - invKmm * DkmmDx{k} * invKmm;
    for u = 1:U
        DknmkmmInvDx = model.Knm{u} * DkmmInvDx + DknmDx{k}{u} * invKmm;
        DknkmInvknDx = DknmkmmInvDx*model.Knm{u}' + KnmInvKmm{u}* DknmDx{k}{u}' ;
        DbDx = diag( (1+1/model.prior.xi)*(DknnDx{k}{u} -DknkmInvknDx)...
            + 2* DknmkmmInvDx * (model.var.Sigma{u} * KnmInvKmm{u}'));
        aDaDx = diag((KnmInvKmm{u}*(model.var.Mu{u}-model.prior.g)+model.prior.g)...
            *(model.var.Mu{u}-model.prior.g)'* DknmkmmInvDx' );
        for i=1:model.Nu(u)
            if (w{u}(i))
                [~, dgg] = queryGz(-meanRate2{u}(i)/2/varRate{u}(i), model.gz, model.dgz);
                df(k) = df(k) - dgg* (-1/2/varRate{u}(i)^2)...
                    *(2*varRate{u}(i)* aDaDx(i)- meanRate2{u}(i)* DbDx(i)) + 1/varRate{u}(i)*DbDx(i);
            end;
        end;
        df(k) = df(k) - v{u}'/model.ratio* (2* aDaDx + DbDx) ...
            + 0.5*trace(invKmm*  (model.var.Sigma{u} + (model.var.Mu{u}-model.global.Mu) * (model.var.Mu{u}-model.global.Mu)' )...
            * invKmm * DkmmDx{k});
    end;
    df(k) = df(k)  - 0.5*(U+1) * trace(invKmm*DkmmDx{k}) + 0.5 * model.prior.xi *...
        trace(invKmm *(model.global.Mu-model.prior.g)*(model.global.Mu-model.prior.g)'*invKmm* DkmmDx{k});
end;
df = -df;
end


