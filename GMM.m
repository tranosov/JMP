function [G]=GMM(pars_,pars,momentest,W,momentall,params)
% set model
% solve model
% create moments to compare to daya

global VERBOSE
if VERBOSE
array2table(pars_, 'RowNames',pars.Properties.RowNames)
end

tic
params(pars.Properties.RowNames,:)= array2table(pars_, 'RowNames',pars.Properties.RowNames);
parsc=calibrate1(params,momentall); %calibrate within
params(parsc.Properties.RowNames,:)=parsc;
[EQS,PARREST] = set_model_estimateR(params,1);
time=toc;

if VERBOSE
toc
end

[p,LA,EXITFLAG,time_]=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0); % just prices! % do not reset params
time=time+time_;


UP=100;
if EXITFLAG==999
    G=999*UP;
else
    FORCEFIT=1; %change THETAs to fit the required marriage behavior
    [moments_,~,time_,EXITFLAG]...
        =moments_withmm(p,LA,params,EQS,PARREST,1,1,FORCEFIT);
    time=time+time_;
    if EXITFLAG==999
        G=999*UP;
    else  
        moments_=moments_(momentest.Properties.RowNames,:);
        G=(table2array(moments_)-table2array(momentest))'*W*(table2array(moments_)-table2array(momentest))*UP;
    end
end

momnames=moments_.Properties.RowNames;
pnames=params.Properties.RowNames;

filename = "./estimation/progress.txt";
io = fopen(filename,'a');
fprintf(io,"GMM function value =   %16.8f\n",G);
fprintf(io," \n");
fprintf(io,"time =   %16.8f\n",time);
fclose(io);


filename = "./estimation/progressmoment.txt";
io = fopen(filename,'a');
fprintf(io,"GMM function value =   %16.8f\n",G);
fprintf(io," \n");
fprintf(io,"time =   %16.8f\n",time);

fprintf(io," data and simulated moments\n");
for jj=1:size(moments_,1)
    fprintf(io,"%s",char(momnames(jj,1)));
    fprintf(io,"%16.8f",table2array(momentest(jj,1)) );  
    fprintf(io,"%16.8f\n",table2array(moments_(jj,1))) ; 
end
fprintf(io," \n");
fprintf(io," parameter value\n");
for jj =1:size(params,1)
    fprintf(io,"%s",char(pnames(jj,1)));
    fprintf(io,"%16.8f\n",table2array(params(jj,1)));
end
fclose(io);





end