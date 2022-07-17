function [VAR,VARt,jac0,pars,parsraw]=SEs(pars_,pars,momentest,W,momentall,forparams,jac0file)
params=forparams(:,'value');
parder=table2array(forparams(forparams.('forses')==1,'deroftransform'));

[moments_,~,~,params_final]=GMMmoments(pars_,pars,momentest,momentall,params,1,0);
m_=moments_;
G=(moments_-table2array(momentest))'*W*(moments_-table2array(momentest))/10^3

[moments_,~,~,params_final2]=GMMmoments(pars_,pars,momentest,momentall,params_final,1,0); % make sure all inputs are proper! (and so is THETAHW)
m2_=moments_; % it is terryfyingly different. working on it
G2=(moments_-table2array(momentest))'*W*(moments_-table2array(momentest))

%mm_=m_-m2_ % it helps a tiny bit to be stricter in labors CONS=100000; but it is super slow and still sucks
 % Only strengthening th eprecision of the multiplier is even better!
parsall=params_final(:,'value');
parsraw=parsall(forparams.('forses')==1,'value');

parsall=transform(parsall); % ALL
pars=parsall(forparams.('forses')==1,'value');
pars_=table2array(pars);


%[moments_,~,~,params_final2]=GMMmoments(pars_,pars,momentest,momentall,parsall,0.5);
%moments_-m_
%table2array(params_final2)-table2array(params_final)


%parsall=params_final2(:,'value');
%parsall(forparams.('forses')==1,'value')
%parsall=transform(parsall); % ALL

%pars=parsall(forparams.('forses')==1,'value')
%pars_=table2array(pars);

moms=@(p) GMMmoments(p,pars,momentest,momentall,parsall,0.5); % lets not recalibrate, except for THETAHW
del_=0.5;
[jac,err,jac0] = jacobianest(moms,pars_,pars.Properties.RowNames,1,jac0file,del_);
G0=jac0' ;
VAR=G0*W*G0'
VAR=VAR\eye(size(VAR,1));
xlswrite(jac0file,VAR,'VAR','A1')


G0=jac0'.*repmat(parder,1,size(jac0',2)) ;
VARt=G0*W*G0'
VARt=VARt\eye(size(VARt,1));
xlswrite(jac0file,VARt,'VARt','A1')


end