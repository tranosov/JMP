function [G]=GMM_noL_MM(pars_,pars,momentest,W,momentall,params,Lpar,Wsq)
pars_=reshape(pars_,max(size(pars_)),1);
% set model
% solve model
% create moments to compare to daya
global VERBOSE GMIN ITER IN %DOWN
global filename1 filename2

if Lpar~=1
    %moments_(end,1)=table2array(params{'LA0',:});% I moved it inside
    %GMMmoments?
%else
    momentest('L',:)=[]; % do not compute L moment
end

    try

        [moments_,time, EXITFLAG,params_]=GMMmoments(pars_,pars,momentest,momentall,params,0.9,0); % still recalibrate (but only sigmam this time)
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

a=1;


if EXITFLAG==999
    params_;
    G=10^(12);
    clear('global', 'IN')
else
    G=(moments_-table2array(momentest))'*W*(moments_-table2array(momentest)); % have to make sure there is W value
    %G=round(G,6); %last instance - get rid of tiny differences!
end


momnames=momentest.Properties.RowNames;
pnames=params_.Properties.RowNames; % params used
pnames2=pars.Properties.RowNames; % estimation inputs

fprintf(" \n");
fprintf("no L:");
fprintf(" \n");
fprintf("GMM function value =   %16.8f\n",G);
fprintf(" \n");
fprintf("time =   %16.8f\n",time);
fprintf(" \n");
fprintf(" iter: %16.8f\n",ITER);

io = fopen(filename1,'a');
if Lpar
    fprintf(io,"GMM function value, L not recomputed=   %16.8f\n",G);
else
    fprintf(io,"GMM function value no L =   %16.8f\n",G);
end

fprintf(io,"time =   %16.8f\n",time);
fprintf(io," iter: %16.8f\n",ITER);
if G<GMIN
        fprintf(" \n");
        fprintf(" FLAG: MINIMUM \n");
        fprintf(io," FLAG: MINIMUM \n");
end

global Wadd GMINadd
if (~isempty(Wadd)) &&  (EXITFLAG<999)
    fprintf(" \n");
    Gadd=(moments_-table2array(momentest))'*Wadd*(moments_-table2array(momentest));
    fprintf(io,"GMM - W=inv(V)   %16.8f\n",Gadd);
    if Gadd<GMINadd
        GMINadd=Gadd;
        fprintf(" \n");
        fprintf(io," FLAG: MINIMUM in W=inv(V) \n");
        fprintf(" \n");
    end

    fprintf(" \n");
end

fclose(io);
ITER=ITER+1;


%filename2 = "./estimation/progressmoment.txt";
io = fopen(filename2,'a');
fprintf(io,"GMM =   %16.8f\n",G);
fprintf(io,"time =   %16.8f\n",time);



if EXITFLAG~=999
    fprintf(" \n");
    fprintf(io," data and simulated moments \n");
    if G<GMIN
        %fprintf(" \n");
        %fprintf(" FLAG: MINIMUM \n");
        fprintf(io," FLAG: MINIMUM \n");
    end
    %weightedm=((moments_-table2array(momentest))'*W).*(moments_-table2array(momentest))';
    %sum(weightedm)
    % shouldn't this pe positive in all elements???
    for jj=1:size(moments_,1)
        fprintf(io,"%s",char(momnames(jj,1)));
        fprintf(io,"%16.8f",table2array(momentest(jj,1)) );  
        fprintf(io,"%16.8f\n",moments_(jj,1)) ; 
    end
    fprintf(io," \n");
    fprintf(io," data and simulated moments: and dif and weighted dif \n");
    weightedm=((moments_-table2array(momentest))'*Wsq); %.*(moments_-table2array(momentest))';
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



if G<GMIN
        GMIN=G;
end

end