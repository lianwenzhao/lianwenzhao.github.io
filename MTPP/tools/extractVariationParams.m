function V = extractVariationParams(model, uOpt)
% index = [];
% for i=1:model.m
%     index = [index; [i:model.m]'+(i-1)*model.m];
% end;
% V = [model.var.Mu{uOpt}; model.var.L{uOpt}(index)];

V =  [model.var.Mu{uOpt}; reshape(model.var.L{uOpt},[model.m^2,1])];
% V = [model.var.Mu{uOpt}; model.var.Sigma{uOpt}(index)];
%[ii,jj] = ind2sub([model.m, model.m],index);
%full(sparse(ii,jj,V(model.m+1:end),model.m,model.m));
