
function locs=sumC(X)
global sizehs

locs=reshape(sum(sum(sum(sum(sum(sum(X,7),6))))),sizehs);
end