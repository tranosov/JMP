%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation. CHECK NEW VERSION');

forparams =readtable('./input/PREP.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

paramsall=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
momentall=formom(:,'value');
momentest=formom(formom.('toestimate')==1,'value');
W=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames))\eye(size(momentest.Properties.RowNames,1));



pars=paramsest;
global RESC
RESC=10^3;
global VERBOSE
VERBOSE=0;

fopt=10^(-6); % no idea
x0=table2array(paramsest);
LB=table2array(forparams(paramsest.Properties.RowNames,'min'));
UB=table2array(forparams(paramsest.Properties.RowNames,'max'));
GG=GMM(x0,pars,momentest,W,momentall,paramsall);


ObjectiveFunction=@(x)GMM(x,pars,momentest,W,momentall,paramsall); % pars has the list!

options = optimoptions(@simulannealbnd,'MaxFunctionEvaluations',10000,'Display','diagnose');
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,x0,LB,UB,options)
output