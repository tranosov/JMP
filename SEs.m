function [SEs]=SEs(pars_,pars,momentest,W,momentall,params)


moms=@(p) GMMmoments(p,pars,momentest,W,momentall,params);
[jac,err] = jacobianest(moms,pars_);
G0=jac ;

%todo: probably rescale parameters so they are all around 1!
% this is insanely slow now. running the routine many times. why?
% I rescaled and still. I mean. it is one derivative at a time, so that is
% not that weird - npar*nmom derivatives. the function is a vector though -
% shouldn't it be just npar number of derivatives?



end