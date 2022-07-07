%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation. CHECK NEW');


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

paramsall=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
momentall=formom(:,'value');
momentest=formom(formom.('toestimate')==1,'value');
W=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames))\eye(size(momentest.Properties.RowNames,1));
%[SEs]=SEs(x0,pars,momentest,W,momentall,paramsall);



pars=paramsest;
global RESC
RESC=10^3;
global VERBOSE
VERBOSE=0;


%fopt=10^(-6); % no idea
x0=table2array(paramsest);
LB=table2array(forparams(paramsest.Properties.RowNames,'min'));
UB=table2array(forparams(paramsest.Properties.RowNames,'max'));
GG=GMM(x0,pars,momentest,W,momentall,paramsall);
GG
%GG=6;


ObjectiveFunction=@(x)GMM(x,pars,momentest,W,momentall,paramsall); % pars has the list!
rng(357);
options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',10000,'Display','diagnose');
% increase temp to have more acceptence
options.InitialTemperature = 100*max(GG,1); %I think should be in scale of objective function (or like jacobian - how params affect obj function)
%temperature = @(optimValues,options) options.InitialTemperature.*(0.99.^optimValues.k); % slow down?
%options.TemperatureFcn=temperature;
options.ReannealInterval=50;

[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,x0,LB,UB,options)
output

% todo: how to make the algorithm look harder?

% maybe lower ReannealInterval a little bit?

% temperature slower?
%temperature = @(optimValues,options) options.InitialTemperature.*(0.99^optimValues.k)
%options.TemperatureFcn=@temperature
