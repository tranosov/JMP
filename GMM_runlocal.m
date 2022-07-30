%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation.\n');


global filename1 filename2
filename1 = "./estimation/progress.txt";
filename2 = "./estimation/progressmoment.txt";
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
%[SEs]=SEs(x0,pars,momentest,W,momentall,paramsall);


momentest_noL=momentest;
momentest_noL('L',:)=[];
W_noL=(DOWN.*table2array(forw( momentest_noL.Properties.RowNames, momentest_noL.Properties.RowNames)))\eye(size(momentest_noL.Properties.RowNames,1));

fprintf('W inverted.\n');


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


global GMIN ITER
GMIN=99999;
ITER=0;


io = fopen(filename1,'a');
fprintf(" \n");
fprintf(io," FLAG: MINIMUM \n");
fclose(io);

% options
options = optimset('Display','iter');
    ObjectiveFunction=@(x)GMM_noL(x,pars,momentest,W_noL,momentall,paramsall); % pars has the list!
    rng(359);
    %options.TolFun=10^(-6);
    %options.TolX=10^(-6);
    %options.MaxFunctionEvaluations

[x,fval,exitFlag,output] = fminsearch(ObjectiveFunction,x0,options);
output

% try fminunc?

io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Routine terminated. \n");
fprintf(io,"%s",output);
fprintf(io," \n");
fclose(io);

% 
