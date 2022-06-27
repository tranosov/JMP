
function locs=sumS(X)
global sizehs

locs=reshape(sum(sum(sum(X,4))),sizehs);
end