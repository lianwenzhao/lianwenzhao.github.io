function model = returnVariationParams(model, V, uOpt)
% index = [];
% for i=1:model.m
%     index = [index; [i:model.m]'+(i-1)*model.m];
% end;
% [ii,jj] = ind2sub([model.m, model.m],index);
% model.var.L{uOpt} = full(sparse(ii,jj,V(model.m+1:end),model.m,model.m));
model.var.L{uOpt} = reshape(V(model.m+1:end), [model.m, model.m]);
model.var.Sigma{uOpt} = model.var.L{uOpt} * model.var.L{uOpt}';
model.var.Mu{uOpt} = V(1:model.m);

