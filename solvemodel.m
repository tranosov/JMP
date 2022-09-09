
function [p,LA,EXITFLAG,time]=solvemodel(params,EQS,PARREST,p0,LA0,reinitialize,psolve,msolve)
% options to kill one or both market clearings?
% output all moments or only

%reinitialize=0.5; % 0 vs 0.5 vs 1 - 0.5 is reset PARREST and EQS with new lambda, but not new inputs  

global RESC
if isempty(RESC)
    RESC=10^6;
end


if reinitialize==0.5
    inputs_=EQS.('inputs');
    [EQS,PARREST]=set_model_estimateR(params,0); % this old version?
    EQS.('inputs')=inputs_;
elseif reinitialize==1
    [EQS,PARREST]=set_model_estimateR(params,1);
end

% assign EQS,PARREST

if psolve==1 && msolve==1
    [x,EXITFLAG,time]=solve([log(p0),LA0*(RESC)],EQS,PARREST);
    p=exp(x(1:end-1));
    LA=x(end)/(RESC);
    
elseif psolve==1 && msolve==0.5 % not lambda, but share of married overall can change
    [x,EXITFLAG,time]=solve_halfmm([log(p0)],EQS,PARREST);
    p=exp(x(1:end));
    LA=LA0;    
elseif psolve && (~msolve)
    [x,EXITFLAG,time]=solvep(log(p0),EQS,PARREST,1);
    p=exp(x);
    LA=LA0;
elseif (~psolve) && msolve
    [x,EXITFLAG,time]=solvem(LA0*(RESC),EQS,PARREST);
    p=p0;
    LA=x/(RESC);
else
    p=p0;
    LA=LA0;
    EXITFLAG=1;
    time=0;
end


end