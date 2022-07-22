%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation.\n');


filename = "./estimation/progressmoment.txt";
io = fopen(filename,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);
filename = "./estimation/progress.txt";
io = fopen(filename,'a');
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
GG=GMM(x0,pars,momentest,W_,momentall,paramsall);
global GMIN ITER
GMIN=GG
ITER=0;
%GG=6;
filename = "./estimation/progress.txt";
io = fopen(filename,'a');
fprintf(" \n");
fprintf(io," FLAG: MINIMUM \n");
fclose(io);

    ObjectiveFunction=@(x)GMM(x,pars,momentest,W_,momentall,paramsall); % pars has the list!
    rng(357);
    options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',10000,'Display','diagnose');
    % increase temp to have more acceptence
    options.InitialTemperature = 1.5*max(GG,1); %no try to bring down, default is 100, I think should be in scale of objective function (or like jacobian - how params affect obj function)
    options.InitialTemperature = 1*max(GG,1); %no try to bring down, default is 100, I think should be in scale of objective function (or like jacobian - how params affect obj function)
    temperature = @(optimValues,options) options.InitialTemperature.*(0.999.^optimValues.k); % slow down?
    options.TemperatureFcn=temperature;
    options.ReannealInterval=1; % brought down A LOT so there is more search? but so far not helping much...
% reannealing is limited IF: temperature already decreased a lot. if
% initial objective is not crazy high - could be more useful?
    
    %options.MaxStallIterations=20; % not sure if this is not just desperate? why would I want to evaluate that many more times around no change?
    %options.FunctionTolerance=10^(-6);

    options.ObjectiveLimit=10^(-6); 
    %options.MaxFunctionEvaluations

[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,x0,LB,UB,options);
output

filename = "./estimation/progress.txt";
io = fopen(filename,'a');
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
