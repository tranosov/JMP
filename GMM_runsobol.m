%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation - sobol points and fminsearch, with L=LA0.\n');


global filename1 filename2
filename1 = "./estimation/progress_sobol.txt";
filename2 = "./estimation/progressmoment_sobol.txt";
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);

io = fopen(filename2,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);

%%
forparams =readtable('./input/PREP.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

fprintf('Inputs loaded.\n');


paramsall=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
momentall=formom(:,'value');
momentest=formom(formom.('toestimate')==1,'value');
global DOWN
DOWN=10^3; % scale down W. it is reallyhigh.

filename = "./estimation/progress.txt";
io = fopen(filename,'a');
fprintf(io," \n");
fprintf(io,"Scaling down weighting matrix by a factor %16.8f\n",DOWN);
fclose(io);


W_=(DOWN.*table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames)))\eye(size(momentest.Properties.RowNames,1));
Wsq=chol(W_);
%[SEs]=SEs(x0,pars,momentest,W,momentall,paramsall);



%momentest_noL=momentest;
%momentest('L',:)=[];
%W_noL=(DOWN.*table2array(forw( momentest_noL.Properties.RowNames, momentest_noL.Properties.RowNames)))\eye(size(momentest_noL.Properties.RowNames,1));

fprintf('W inverted.\n');

%paramsest('LA0',:)=[];
pars=paramsest;
global RESC
RESC=10^3;
global VERBOSE
VERBOSE=0;


%fopt=10^(-6); % no idea
x0=table2array(paramsest);
LB=table2array(forparams(paramsest.Properties.RowNames,'min'));
UB=table2array(forparams(paramsest.Properties.RowNames,'max'));
%GGnoL=GMM_noL(x0,pars,momentest,W_noL,momentall,paramsall);


% run a short version from my x0
global GMIN ITER
GMIN=99999;
ITER=0;


io = fopen(filename1,'a');
fprintf(" \n");
fprintf(io," FLAG: MINIMUM \n");
fprintf(io,"GMM =   %16.8f\n",GMIN);
fclose(io);
% options
options = optimset('Display','iter');
    ObjectiveFunction=@(x)GMM_noL(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
    rng(354);
    options.TolFun=10^(-1); % lenient!
    options.TolX=10^(0); % hopefully  - should not matetr
    options.MaxFunctionEvaluations=500; % pretty low!

    %Unlike other solvers, fminsearch stops when it satisfies both TolFun and TolX.

[x,fval,exitFlag,output] = fminsearch(ObjectiveFunction,x0,options);
output


io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Routine terminated. \n");
fprintf(io,"%s",output);
fprintf(io," \n");
fclose(io);




NP=size(table2array(pars),1);
sob = sobolset(NP,'Skip',1e3,'Leap',1e2);
sob= scramble(sob,'MatousekAffineOwen');
SOBN=10;
sob_ = sob(1:SOBN,:);
for sobn=1:SOBN
    io = fopen(filename1,'a');
    fprintf(io," \n");
    fprintf(io,"New input0. \n");
    fclose(io);

    io = fopen(filename2,'a');
    fprintf(io," \n");
    fprintf(io,"New input0. \n");
    fclose(io);

    x00=UB.*sob_(1,:)'+LB.*(1-sob_(1,:))';
    % options
    options = optimset('Display','iter');
        ObjectiveFunction=@(x)GMM_noL(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
        rng(354);
        options.TolFun=10^(-1); % lenient!
        options.TolX=10^(0); % hopefully  - should not matetr
        options.MaxFunctionEvaluations=2000; % pretty low!

        %Unlike other solvers, fminsearch stops when it satisfies both TolFun and TolX.
        
    [x,fval,exitFlag,output] = fminsearch(ObjectiveFunction,x00,options);
    output


    io = fopen(filename1,'a');
    fprintf(io," \n");
    fprintf(io,"Routine terminated. \n");
    fprintf(io,"%s",output);
    fprintf(io," \n");
    fclose(io);

    
end




% try fminunc?

% 1191 iterations and still did not fucking converge. still is improving.
% that is wild.
