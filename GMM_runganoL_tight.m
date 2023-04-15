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
filename1 = "./estimation/progress_ganoL_econstat6_tight_nW.txt";
filename2 = "./estimation/progressmoment_ganoL_econstat6_tight_nW.txt";
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine.  \n");
fclose(io);

io = fopen(filename2,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine.\n");
fclose(io);

%%
forparams =readtable('./input/PREP_tight.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

fprintf('Inputs loaded.\n');

paramsall=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
momentall=formom(:,'value');
momentest=formom(formom.('toestimate')==1,'value');
blow=formom(formom.('toestimate')==1,'BLOW');

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
    blow('clmm',:)={1};
    Wall(end+1,end+1)=DOWN;

Wsq_all=chol(Wall);
Wall2=Wall.*(ones(size(Wall))+diag(table2array(blow))-(diag(ones(size(Wall,1),1))) );
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
x0=table2array(paramsest);
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
GMINadd=999; %G_W;

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

    ObjectiveFunction=@(x)GMM_noL_MM(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
    rng(300);
    options = optimoptions(@ga,'Display','iter');
    % increase temp to have more acceptence
    
intcon = [];
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
