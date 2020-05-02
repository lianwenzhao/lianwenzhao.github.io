function [rate] = feature2rate(feature, ratePar, u, fnProvided)
if (fnProvided==0)
    if strncmp(ratePar.opt, 'pseudo',2)
        if iscell(ratePar.cov{2}{1});
            Knm = feval( ratePar.cov{2}{1}{:},ratePar.GP.logtheta(1:end-1), feature', ratePar.Xu);
        else     Knm = feval( ratePar.cov{2}{1},ratePar.GP.logtheta(1:end-1), feature', ratePar.Xu);
        end;
        rate = (ratePar.prior.g + Knm*ratePar.KinvF{u}).^2;
%         rate = rate +0.5* ...
%             diag(ratePar.diagK - Knm* ratePar.KmmInv * Knm');
%         rate = (ratePar.prior.g + Knm*ratePar.KinvF{u}).^2 +0.5* ...
%             diag(ratePar.diagK - Knm* ratePar.KmmInv * Knm' + Knm*ratePar.KmmInv*ratePar.SigmaInv{u}*ratePar.KmmInv*Knm');
        rate = rate/ratePar.ratio;
        if rate>ratePar.maxRate
            rate = ratePar.maxRate;
        elseif rate<ratePar.minRate
            rate = ratePar.minRate;
        end;
    end;
else
    if strncmp(ratePar.opt,'sigmoid',1)
        rate = feature'*ratePar.weights(:,u);
        rate = ratePar.minRate + ratePar.scaling * 1/(1+exp(-rate));
    elseif strncmp(ratePar.opt, 'pseudo',2)
        if iscell(ratePar.cov{2}{1});
            Knm = feval( ratePar.cov{2}{1}{:},ratePar.GP.logtheta(1:end-1), feature', ratePar.Xu);
        else     Knm = feval( ratePar.cov{2}{1},ratePar.GP.logtheta(1:end-1), feature', ratePar.Xu);
        end;
        rate = (ratePar.fnTable{u}(num2str(feature'))).^2/ratePar.ratio;
        %         rate = exp(Knm*ratePar.KinvF{u} + ratePar.randTable(num2str(feature')) * sqrt(ratePar.diagK - Knm* ratePar.KmmInv * Knm'));% randn()
        if rate>ratePar.maxRate
            rate = ratePar.maxRate;
        elseif rate<ratePar.minRate
            rate = ratePar.minRate;
        end;
        %     rate = max(rate, ratePar.minRate);
        %     rate = min(rate, ratePar.maxRate);
    else
        rate = feature'*ratePar.weights(:,u);
        rate = max(rate, ratePar.minRate);
        rate = min(rate, ratePar.maxRate);
    end;
end;
