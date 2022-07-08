function [G]=GMM(pars_,pars,momentest,W,momentall,params)
% set model
% solve model
% create moments to compare to daya
global VERBOSE %DOWN

[moments_,time, EXITFLAG,params_]=GMMmoments(pars_,pars,momentest,W,momentall,params);


if EXITFLAG==999
    G=10^(23);
else
    G=(moments_-table2array(momentest))'*W*(moments_-table2array(momentest));
end
momnames=momentest.Properties.RowNames;
pnames=params_.Properties.RowNames;



fprintf("GMM function value =   %16.8f\n",G);
fprintf(" \n");
fprintf("time =   %16.8f\n",time);

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


if EXITFLAG~=999
    fprintf(io," data and simulated moments, and dif, absdif, weighted dif \n");
    weightedm=((moments_-table2array(momentest))'*W).*(moments_-table2array(momentest))';
    for jj=1:size(moments_,1)
        fprintf(io,"%s",char(momnames(jj,1)));
        fprintf(io,"%16.8f",table2array(momentest(jj,1)) );  
        fprintf(io,"%16.8f\n",moments_(jj,1)) ; 
    end
    fprintf(io," data and simulated moments: and dif and weighted dif \n");
    weightedm=((moments_-table2array(momentest))'*W).*(moments_-table2array(momentest))';
    for jj=1:size(moments_,1)
        fprintf(io,"%s",char(momnames(jj,1)));
        fprintf(io,"%16.8f",table2array(momentest(jj,1)) -moments_(jj,1)) ; 
        %fprintf(io,"%16.8f ",abs(table2array(momentest(jj,1)) -moments_(jj,1))) ; 
        fprintf(io,"%16.8f\n",weightedm(1,jj)) ; 
    end
else
    fprintf(io," data and simulated moments skipped.\n");
    fprintf(io," solvemodel failed.\n");
end

fprintf(io," \n");
fprintf(io," parameter value\n");
for jj =1:size(params,1)
    fprintf(io,"%s",char(pnames(jj,1)));
    fprintf(io,"%16.8f\n",table2array(params_(jj,1)));
end
fclose(io);





end