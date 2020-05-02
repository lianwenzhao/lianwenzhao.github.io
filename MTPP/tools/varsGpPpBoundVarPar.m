function [f, df] = varsGpPpBoundVarPar(W, model, uOpt, flag)
%%%%compute f and df of lower bound w.r.t. variational parameters

% PRELIMINARIES
if nargin == 2
    flag = 1;
end

% place W to the model structure
model = returnVariationParams(model, W, uOpt);

% U = length(model.events);

% COVARIANCE FUNCTION QUANTITIES needed in the sparse method: K_mm, Knm, diag(Knn)
%
% inducing variable parameters are pseudo-inputs
if strcmp(model.indType, 'pseudoIns')
    model.Kmm = feval(model.cov{:}, model.GP.logtheta,model.Xu);
    for u = uOpt:uOpt% 1:U
        if iscell(model.cov{2}{1});
            model.Knm{u} = feval(model.cov{2}{1}{:},model.GP.logtheta(1:end-1),model.events{u}.feature,model.Xu);
        else
            model.Knm{u} = feval(model.cov{2}{1},model.GP.logtheta(1:end-1),model.events{u}.feature,model.Xu);
        end;
        %         model.Knn{u} = feval(cov{:},model.GP.logtheta,model.events{u}.feature);
    end;
    % inducing variables are linear combinations of training latent variables
elseif strcmp(model.indType, 'weights')
    %%%%%%%%%not implemented yet!!!!
    % only inducing variables are optimized
    % this allows kernel quantities to be precomputed once
    if flag==2 & isfield(model, 'KmmSq')
        model.Kmm = zeros(model.m, model.m);
        model.Knm = zeros(model.n, model.m);
        cnt = 0;
        for j=1:model.nIndParams
            W = sparse(diag(model.Xu(:,j)));
            model.Kmm = model.Kmm + W*(model.KmmSq(:,:,j)*W);
            model.Knm = model.Knm + model.KnmAll(:,:,j)*W;
            % cross-terms for K_mm
            for k=j+1:model.nIndParams
                cnt = cnt + 1;
                G = sparse(diag( model.Xu(:,j) ))*( model.KmmCr(:,:,cnt)*sparse(diag(model.Xu(:,k))) );
                model.Kmm = model.Kmm + (G + G');
            end
        end
    else
        [model.Kmm model.KmmSq model.KmmCr] = kernelWeights(model);
        [model.Knm model.KnmAll] = kernelWeights(model, model.X);
    end
end

model.diagKnn = sum(exp(2*model.GP.logtheta(end-1:end)));%%%simply the sigma_f^2
if model.GP.constDiag == 1
    for u=uOpt:uOpt%1:U
        model.TrKnn{u} = model.Nu(u)*model.diagKnn;
    end;
else
    error('Knn not properly defined! ');
end


% upper triangular Cholesky decomposition
% (we add jitter to Kmm which implies jitter inducing variables
% however the matrix that is stored is jitter free.
% The jitter-free matrix is used to compute more precise derivatives; see
% documentation)
for u = uOpt:uOpt%1:U
    w{u} = model.events{u}.label;
    v{u} = model.events{u}.duration;
end;
% invKmm = (model.Kmm+model.jitter*eye(model.m))\eye(model.m);
invKmm = (model.Kmm)\eye(model.m);
f = 0;
for u = uOpt:uOpt%1:U
    Dmean = - model.global.invSigma * (model.var.Mu{u}-model.global.Mu);
    Linv = model.var.L{u}\eye(model.m);
    Dvar = -0.5* model.global.invSigma + 0.5* (Linv'*Linv);
%     Dvar = -0.5* model.global.invSigma + 0.5* model.var.Sigma{u}\eye(model.m);
    KnmInvKmm{u} = model.Knm{u}*invKmm;
    meanRate = KnmInvKmm{u}*(model.var.Mu{u}-model.prior.g) + model.prior.g;
    meanRate2 = meanRate.^2;
    varRate = (1+1/model.prior.xi)*(model.diagKnn - diag(KnmInvKmm{u}*model.Knm{u}'))...
        + diag(KnmInvKmm{u}*model.var.Sigma{u}*KnmInvKmm{u}');
    for i=1:model.Nu(u)
        if w{u}(i)
           [gg, dgg] = queryGz(-meanRate2(i)/2/varRate(i), model.gz, model.dgz);
            f = f -  gg+ log(0.5*varRate(i));
            Dmean = Dmean + dgg/varRate(i) * meanRate(i)* KnmInvKmm{u}(i,:)';
            Dvar = Dvar + (-dgg*meanRate2(i)/2/varRate(i)^2+ 1/varRate(i)  ) * (KnmInvKmm{u}(i,:)' *KnmInvKmm{u}(i,:));
        end;
        Dvar = Dvar - v{u}(i)/model.ratio*(KnmInvKmm{u}(i,:)' *KnmInvKmm{u}(i,:));
    end;
%     disp([f, - v{u}'*(meanRate2+varRate), 0.5*logdet(model.var.Sigma{u})...
%         - 0.5* trace(model.global.invSigma * (model.var.Sigma{u} + ...
%         (model.var.Mu{u} - model.global.Mu)*(model.var.Mu{u} - model.global.Mu)' )), ...
%         sum([f, - v{u}'*(meanRate2+varRate), 0.5*logdet(model.var.Sigma{u})...
%         - 0.5* trace(model.global.invSigma * (model.var.Sigma{u} + ...
%         (model.var.Mu{u} - model.global.Mu)*(model.var.Mu{u} - model.global.Mu)' )) ])])
%     disp([model.var.Mu{u} Dmean+model.global.invSigma * (model.var.Mu{u}-model.global.Mu), - 2* (v{u}'*diag(meanRate)*KnmInvKmm{u})',...
%         model.global.invSigma * (model.var.Mu{u}-model.global.Mu),...
%         sum([Dmean+model.global.invSigma * (model.var.Mu{u}-model.global.Mu), - 2* (v{u}'*diag(meanRate)*KnmInvKmm{u})',...
%         model.global.invSigma * (model.var.Mu{u}-model.global.Mu)],2)])
    f = f - v{u}'/model.ratio*(meanRate2+varRate);
    Dmean = Dmean - 2* ((v{u}/model.ratio.*meanRate)'*KnmInvKmm{u})';
    f = f  + 0.5*logdet(model.var.Sigma{u})...
        - 0.5* trace(model.global.invSigma * (model.var.Sigma{u} + ...
        (model.var.Mu{u} - model.global.Mu)*(model.var.Mu{u} - model.global.Mu)' ));
end;

f = -f;
updated.var.Mu{1} = Dmean;
% updated.var.Sigma{1} = 2* Dvar*model.var.L{uOpt};
updated.var.L{1} = 2* Dvar*model.var.L{uOpt};
updated.m = model.m;
df = extractVariationParams(updated, 1);
df = -df;


