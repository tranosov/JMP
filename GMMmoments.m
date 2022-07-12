function [moments_,time,EXITFLAG,params_final]=GMMmoments(pars_,pars,momentest,momentall,params,recalibrate)

%recalibrate: 1 - all within estimation recalibration. 0.6 - only THETAHW THETA, 0.5 - only THETAHW

% set model
% solve model
% create moments to compare to daya

global VERBOSE

tic
params(pars.Properties.RowNames,:)= array2table(pars_, 'RowNames',pars.Properties.RowNames);
params=untransform(params);

%params(pars.Properties.RowNames,:) % already untransformed
if recalibrate==1
    parsc=calibrate1(params,momentall); %calibrate within
    params(parsc.Properties.RowNames,:)=parsc;
end



[EQS,PARREST] = set_model_estimateR(params,1);
params=PARREST.('params');
if VERBOSE
    params %(pars.Properties.RowNames,:) % show untransofmed version!
%array2table(pars_, 'RowNames',pars.Properties.RowNames)
end

time=toc;

if VERBOSE
toc
end

if recalibrate>=0.6 % THETAHW and THETA will be recalibrated, so do not resolve lambda 
[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0); % just prices! % do not reset params

elseif recalibrate==0.5 % only THETAHW should be recalibrated
[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0.5);

else
    % so I have to solve with marriage market to have a new lambda - if I
    % wanted comp stats wrt THETAHW
[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,1);
end

time=time+time_;



if EXITFLAG==999
    moments_=999;
    params_final=params;
    return 
else
    FORCEFIT=(recalibrate>=0.6)+0.5*(recalibrate==0.5); %change THETAs to fit the required marriage behavior?
    [moments_,~,time_,EXITFLAG,params_final]...
        =moments_withmm(p,LA,params,EQS,PARREST,1,1,FORCEFIT);
    time=time+time_;
    if EXITFLAG==999
        return
    else  
        moments_=moments_(momentest.Properties.RowNames,:);
    end


end
moments_=table2array(moments_);
moments_(isnan(moments_))=999*10^3; % this is pretty arbitrary. is primarily because sometimes the model has now pnw. And I guess I want them?


end