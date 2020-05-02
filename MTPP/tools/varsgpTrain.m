function [model, margLogL] = varsgpTrain(model, options)
% Inputs:
%         * model: Structure containing the GP model; see varsgpCreate.
%         * options: Several options used in the optimization. options(1)
%           is the maximum number of the objective function evaluations.
%           options(2) is the type of optimization: at the moment only scaled
%           conjugate gradients is used and particularly the minimize.m function of
%           Carl Rasmussen. options(3) = 1,2,or 3; when 1 we optimize over both
%           model (kernel and likelihood) hyperoarameters and inducing variables parameters,
%           when 2 we optimize only over inducing variables parameters, when 3 we optimize
%           over only model hyperparameters.
% Outputs:
%         * model: the model structure with the optimized parameters
%         * margLogL:  the final value of the variational lower bound on the
%           log marginal likelihood.

FuncEval = -options(1);

U = length(model.Nu);
[AA, bb, AAe, bbe] = cholConstraint(model.m);
fEM = [];
fM = 0;
fE1 = 0;
for u = 1:U
    V = extractVariationParams(model, u);
    myObj = @(VV)varsGpPpBoundVarPar(VV, model, u, options(3));
    optionsOpt = optimset('GradObj','on','GradConstr','on','Display','off','Algorithm', 'active-set');
    [V, ftmp] = fmincon(myObj, V, AA, bb, AAe, bbe,[],[],[], optionsOpt);
    model = returnVariationParams(model, V, u);
    fE1 = fE1 - ftmp(end);
end;

fPre = -1e10;
for it = 1:options(4)
    if strcmp(model.objectFunc, 'var')
        % extract the parameters to be optimized
        if (options(3)~=0)
            W = extractOptimizedParams(model, options(3));
            [W, fM] = minimize(W, 'varsGpPpBoundHyper', FuncEval, model, options(3));
            model = returnOptimizedParams(model, W, options(3));
        end;
        [model, fE2]= varsGpPpGlobalPar(model);
        fE1 = 0;
        for u = 1:U
            V = extractVariationParams(model, u);
            myObj = @(VV)varsGpPpBoundVarPar(VV, model, u, options(3));
            optionsOpt = optimset('GradObj','on','GradConstr','on','Display','off','Algorithm', 'active-set');
            [V, ftmp] = fmincon(myObj, V, AA, bb, [],[],[],[],[], optionsOpt);
            model = returnVariationParams(model, V, u);
            fE1 = fE1 - ftmp(end);
        end;
        fEM = [ fEM, fE1 - fE2 - fM(end)];
    end
    disp(['Iterartion ' num2str(it) ' E ' num2str(fEM(end))]);%' M ' num2str(fM(end))
%     save(['./results/model' num2str(it)],'model');
    if abs((fEM(end)-fPre)/fPre)<1e-6;
        break;
    end;
    fPre = fEM(end);
end;
margLogL = fEM;
