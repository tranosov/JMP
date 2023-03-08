%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation - sm, with L=LA0.\n');


global filename1 filename2
filename1 = "./estimation/progress_snnoL_econstat5.txt";
filename2 = "./estimation/progressmoment_smnoL_econstat5.txt";
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine.  \n");
fclose(io);

io = fopen(filename2,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine.\n");
fclose(io);

%%
forparams =readtable('./input/PREP.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
forparams_x0 =readtable('./input/PREP.xlsx','Sheet','PARSTEST','ReadVariableNames', true,'ReadRowNames',true);
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

fprintf('Inputs loaded.\n');


paramsall=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
paramsest_x0=forparams_x0(forparams.('toestimateR')==1,'value');
momentall=formom(:,'value');
momentest=formom(formom.('toestimate')==1,'value');

global DOWN
DOWN=10^3;
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Scaling down weighting matrix by a factor %16.8f\n",DOWN);
fclose(io);

Wall=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames));
Wall=(DOWN.*Wall)\eye(size(momentest.Properties.RowNames,1));
% ADD MM CLEARING CONDITION
    momentest('clmm',:)={0}; 
    Wall(end+1,end+1)=DOWN; 
 
Wsq_all=chol(Wall);
select=momentest;
select.('blowW')=ones(size(Wall,1),1);
select('L','blowW')={1};
select('scommiles','blowW')={10^2};
%select('swcommiles_difw','blowW')={10^3};
select('shcommiles_dif','blowW')={10^3};
select('hdj','blowW')={10^2};
select('sdj_dif','blowW')={10^2};
select('wagegap_hw_withn','blowW')={10^2};
select('betahrs_w','blowW')={10^2};
select('p_gradient_simple','blowW')={10^2};
%select('betalwg_w','blowW')={1};
Wall2=Wall.*diag(select.('blowW'));
Wsq_all2=chol(Wall2); %Wsq_all.*diag(select.('blowW'));


%W_=(DOWN.*table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames)))\eye(size(momentest.Properties.RowNames,1));
%Wsq=chol(W_);

%momentest_noL=momentest;
%momentest('L',:)=[];
%W_noL=(DOWN.*table2array(forw( momentest_noL.Properties.RowNames, momentest_noL.Properties.RowNames)))\eye(size(momentest_noL.Properties.RowNames,1));

fprintf('W inverted.\n');

%paramsest('LA0',:)=[];
pars=paramsest;
global RESC
RESC=10^3;

%fopt=10^(-6); % no idea

x0=table2array(paramsest_x0);
LB=table2array(forparams(paramsest.Properties.RowNames,'min'));
UB=table2array(forparams(paramsest.Properties.RowNames,'max'));

%G_Wdiag=GMM_noL(x0,pars,momentest,Wdiag,momentall,paramsall,1,Wsq_diag);
%G_W=GMM_noL_MM(x0,pars,momentest,Wall,momentall,paramsall,1,Wsq_all);

%G_Wdiag2=GMM_noL(x0,pars,momentest,Wdiag2,momentall,paramsall,1,Wsq_diag2);
%G_W2=GMM_noL_MM(x0,pars,momentest,Wall2,momentall,paramsall,1,Wsq_all2);

W_=Wall2;
Wsq=Wsq_all2;

global Wadd GMINadd
Wadd=Wall;
GMINadd=999999999999999; %G_W;

global VERBOSE
VERBOSE=0;


% run a short version from my x0
global GMIN ITER
GMIN=99999;
ITER=0;


%%

io = fopen(filename1,'a');
fprintf(" \n");
fprintf(io," FLAG: MINIMUM \n");
fclose(io);

G_min=GMM_noL_MM(LB,pars,momentest,W_,momentall,paramsall,1,Wsq);
G_max=GMM_noL_MM(UB,pars,momentest,W_,momentall,paramsall,1,Wsq);


    ObjectiveFunction=@(x)GMM_noL_MM(x,pars,momentest,W_,momentall,paramsall,1,Wsq);
    GG=ObjectiveFunction(x0);
    options.InitialTemperature = 1*max(GG,1); %no try to bring down, default is 100, I think should be in scale of objective function (or like jacobian - how params affect obj function)
    temperature = @(optimValues,options) options.InitialTemperature.*(0.999.^optimValues.k); % slow down?
    options.TemperatureFcn=temperature;
    %options.ReannealInterval=1; % brought down A LOT so there is more search? but so far not helping much...
% reannealing is limited IF: temperature already decreased a lot. if
% initial objective is not crazy high - could be more useful?
    
    options.MaxStallIterations=150; 
    % not sure if this is not just desperate? why would I want to evaluate that many more times around no change?
    options.FunctionTolerance=10^(-3);
    %options.HybridFcn='fminsearch';

    %options.ObjectiveLimit=10^(-6); 
    %options.MaxFunctionEvaluations

[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,x0,LB,UB,options);
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
