function [moments_,time,EXITFLAG]=GMMmoments(pars_,pars,momentest,W,momentall,params)
% set model
% solve model
% create moments to compare to daya

global VERBOSE
if VERBOSE
array2table(pars_, 'RowNames',pars.Properties.RowNames)
end

tic
params(pars.Properties.RowNames,:)= array2table(pars_, 'RowNames',pars.Properties.RowNames);
params=untransform(params);

%params(pars.Properties.RowNames,:) % already untransformed

parsc=calibrate1(params,momentall); %calibrate within
params(parsc.Properties.RowNames,:)=parsc;
[EQS,PARREST] = set_model_estimateR(params,1);
time=toc;

if VERBOSE
toc
end

[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0); % just prices! % do not reset params
time=time+time_;



if EXITFLAG==999
    moments_=999;
    return 
else
    FORCEFIT=1; %change THETAs to fit the required marriage behavior
    [moments_,~,time_,EXITFLAG]...
        =moments_withmm(p,LA,params,EQS,PARREST,1,1,FORCEFIT);
    time=time+time_;
    if EXITFLAG==999
        return
    else  
        moments_=moments_(momentest.Properties.RowNames,:);
    end


end
moments_=table2array(moments_);

end