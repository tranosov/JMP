function [moments_,time,EXITFLAG,params_final]=GMMmoments(pars_,pars,momentest,momentall,params,recalibrate,withmm)

%recalibrate: 1 - all within estimation recalibration. 0.6 - only THETAHW THETA, 0.5 - only THETAHW, 0.1 - only THETA
%withmm: 0 is for not computing lambda moments numerically.
% set model
% solve model
% create moments to compare to daya

global VERBOSE show

tic
params(pars.Properties.RowNames,:)= array2table( reshape(pars_,max(size(pars_)),1), 'RowNames',pars.Properties.RowNames);

if recalibrate>=0.9
    params=untransform(params);
end

%params(pars.Properties.RowNames,:) % already untransformed
if recalibrate>=0.9
    parsc=calibrate1(params,momentall); %calibrate within
    params(parsc.Properties.RowNames,:)=parsc;
end



[EQS,PARREST] = set_model_estimateR(params,1);
params=PARREST.('params');


time=toc;

if VERBOSE
toc
end

if recalibrate>=0.6 | (recalibrate==0.1)  % THETAHW and THETA will be recalibrated oe other reason, so do not resolve lambda 
[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0); % just prices! % do not reset params

elseif recalibrate==0.5 % only THETAHW should be recalibrated
[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0.5);


else
    % so I have to solve with marriage market to have a new lambda - if I
    % wanted comp stats wrt THETAHW
[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,1);
end

time=time+time_;


%if withmm==0
%    if sum(strcmp('L',momentest.Properties.RowNames))>0
%        momentest('L',:)=[];
%    end
%end

if EXITFLAG==999
    moments_=999;
    params_final=params;
    return 
else

    FORCEFIT=(recalibrate>=0.99)+ 0.6*(recalibrate==0.6)+0.5*(recalibrate==0.5)+ 0.1*(recalibrate==0.1)+ 0.9*(recalibrate==0.9); %change THETAs to fit the required marriage behavior?
    [moments_,~,time_,EXITFLAG,params_final]...
        =moments_withmm(p,LA,params,EQS,PARREST,withmm,1,FORCEFIT);
    time=time+time_;
    if EXITFLAG==999
        return
    else  
        moments_=moments_(momentest.Properties.RowNames,:);
    end


end


moments_=table2array(moments_);
moments_(isnan(moments_))=10^12; % this is pretty arbitrary. is primarily because sometimes the model has now pnw. And I guess I want them?


end