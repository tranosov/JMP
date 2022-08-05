%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation - ga, with L=LA0.\n');


global filename1 filename2
filename1 = "./estimation/progress_ganoL.txt";
filename2 = "./estimation/progressmoment_ganoL.txt";
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. new W! \n");
fclose(io);

io = fopen(filename2,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. new W!\n");
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

%filename = "./estimation/progress.txt";
io = fopen(filename1,'a');
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


global GMIN ITER
GMIN=99999;
ITER=0;


io = fopen(filename1,'a');
fprintf(" \n");
fprintf(io," FLAG: MINIMUM \n");
fclose(io);


    ObjectiveFunction=@(x)GMM_noL(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
    rng(300);
    options = optimoptions(@ga,'Display','iter');
    % increase temp to have more acceptence
    
intcon = 1;
nonlcon = [];
A = [];
b = [];
Aeq = [];
beq = [];   
[x,fval,exitFlag,output] = ga(ObjectiveFunction,size(x0,1),A,b,Aeq,beq,LB,UB,nonlcon,intcon,options)
%simulannealbnd(ObjectiveFunction,x0,LB,UB,options);
output


io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Routine terminated. \n");
fprintf(io,"%s",output);
fprintf(io," \n");
fclose(io);

% todo: how to make the algorithm look harder?

% maybe lower ReannealInterval a little bit?

% temperature slower?
%temperature = @(optimValues,options) options.InitialTemperature.*(0.99^optimValues.k)
%options.TemperatureFcn=@temperature
