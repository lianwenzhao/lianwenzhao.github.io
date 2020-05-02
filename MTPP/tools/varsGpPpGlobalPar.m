function [model, f] = varsGpPpGlobalPar(model)
U = length(model.Nu);
model.global.xi = model.prior.xi + U;
model.global.nu = model.prior.nu + U;
model.global.Mu = 1/model.global.xi * sum(cell2mat(model.var.Mu),2) ...
    + model.prior.xi/model.global.xi * model.prior.g;
model.Kmm = feval(model.cov{:}, model.GP.logtheta,model.Xu);
model.global.Sigma = model.Kmm;
model.global.invSigma = model.global.Sigma\eye(model.m);
f = - 0.5 * trace(model.global.invSigma  * ...
    ( ( model.prior.xi * (model.global.Mu-model.prior.g) * (model.global.Mu-model.prior.g)' )));
f = -f;


