function [G]=GMM(pars_,pars,momentest,W,momentall,params)
% set model
% solve model
% create moments to compare to daya
global VERBOSE GMIN ITER IN %DOWN
global filename1 filename2


    try
       %Error-maker
       [moments_,time, EXITFLAG,params_]=GMMmoments(pars_,pars,momentest,momentall,params,1,1);
    catch e %e is an MException struct
        fprintf('The identifier was:\n%s',e.identifier);
        fprintf('There was an error! The message was:\n%s',e.message);
            io = fopen(filename1,'a');
            fprintf(io,'The identifier was:\n%s',e.identifier);
            fprintf(io,'There was an error! The message was:\n%s',e.message);
            % more error handling...
            fclose(io);
        EXITFLAG=999;
        time=0;
        params_=params;
    end



if EXITFLAG==999
    G=10^(6);
    clear('global', 'IN')
elseif EXITFLAG==991
    G=10^(3) + (moments_-table2array(momentest))'*W*(moments_-table2array(momentest));
    clear('global', 'IN')
else
    G=(moments_-table2array(momentest))'*W*(moments_-table2array(momentest));
    G=round(G,6); %last instance - get rid of tiny differences!
end
momnames=momentest.Properties.RowNames;
pnames=params_.Properties.RowNames; % params used
pnames2=pars.Properties.RowNames; % estimation inputs


fprintf(" \n");
fprintf("GMM function value with L =   %16.8f\n",G);
fprintf(" \n");
fprintf("time =   %16.8f\n",time);
fprintf(" \n");
fprintf(" iter: %16.8f\n",ITER);

%filename1 = "./estimation/progress.txt";
io = fopen(filename1,'a');
fprintf(io,"GMM function value =   %16.8f\n",G);
fprintf(io,"time =   %16.8f\n",time);
fprintf(io," iter: %16.8f\n",ITER);
if G<GMIN
        GMIN=G;
        fprintf(" \n");
        fprintf(" FLAG: MINIMUM \n");
        fprintf(io," FLAG: MINIMUM \n");
end

fclose(io);
ITER=ITER+1;


%filename2 = "./estimation/progressmoment.txt";
io = fopen(filename2,'a');
fprintf(io,"GMM function value =   %16.8f\n",G);
fprintf(io,"time =   %16.8f\n",time);



if EXITFLAG~=999
    fprintf(" \n");
    fprintf(io," data and simulated moments \n");
    if G<GMIN
        GMIN=G;
        %fprintf(" \n");
        %fprintf(" FLAG: MINIMUM \n");
        fprintf(io," FLAG: MINIMUM \n");
    end
    weightedm=((moments_-table2array(momentest))'*W).*(moments_-table2array(momentest))';
    sum(weightedm)
    % shouldn't this pe positive in all elements???
    for jj=1:size(moments_,1)
        fprintf(io,"%s",char(momnames(jj,1)));
        fprintf(io,"%16.8f",table2array(momentest(jj,1)) );  
        fprintf(io,"%16.8f\n",moments_(jj,1)) ; 
    end
    fprintf(io," \n");
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
    fprintf(io,"%16.10f\n",table2array(params_(jj,1)));
end
fprintf(io," \n");
fprintf(io," parameter value - raw inputs\n");
for jj =1:size(pars_,1)
    fprintf(io,"%s",char(pnames2(jj,1)));
    fprintf(io,"  ");
    fprintf(io,"%16.15f\n",pars_(jj,1));
end

fprintf(io," \n");
fclose(io);





end