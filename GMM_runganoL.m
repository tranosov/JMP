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
filename1 = "./estimation/progress_ganoL_econstat5.txt";
filename2 = "./estimation/progressmoment_ganoL_econstat5.txt";
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
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

fprintf('Inputs loaded.\n');


paramsall=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
momentall=formom(:,'value');
momentest=formom(formom.('toestimate')==1,'value');
global DOWN
DOWN=1; %10^3; % scale down W. it is reallyhigh. now that I simplify W - get rid of


io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Scaling down weighting matrix by a factor %16.8f\n",DOWN);
fclose(io);


DOWN=10^3;
W__=diag(diag(table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames))));
Wall=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames));
Wall=(DOWN.*Wall)\eye(size(momentest.Properties.RowNames,1));
Wsq_all=chol(Wall);

%[SEs]=SEs(x0,pars,momentest,W,momentall,paramsall);

DOWN=10^3;
Wall=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames));
Wall=(DOWN.*Wall)\eye(size(momentest.Properties.RowNames,1));
Wsq_all=chol(Wall);
select=momentest;
select.('blowW')=ones(size(W__,1),1);
select('L','blowW')={1};
select('scommiles','blowW')={10^4};
%select('swcommiles_difw','blowW')={10^3};
select('shcommiles_dif','blowW')={10^4};
select('hdj','blowW')={10^4};
select('sdj_dif','blowW')={10^4};
select('wlfp_dif','blowW')={10^4};
select('wagegap_hw_withn','blowW')={10^4};
select('betahrs_w','blowW')={10^4};
select('betalwg_w','blowW')={10^4};
select('p_gradient_simple','blowW')={10^4};
%select('betalwg_w','blowW')={1};
Wall2=Wall+diag(select.('blowW')).*diag(Wall); % NOTICE - THERE WAS AN ERROR HERE. NOT ELEMENT BY ELEMNT, HAS TO BE MATRIX MULTIPLICATION. FUCK
Wsq_all2=chol(Wall2);

% does not work - I was doing stupid stuff, but what I had in mind is not
% positive definite it seems.  - lets try an addition????? maybe not enough
% but can try I guess for now?

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


G_W=GMM_noL(x0,pars,momentest,Wall,momentall,paramsall,1,Wsq_all);

%G_Wdiag2=GMM_noL(x0,pars,momentest,Wdiag2,momentall,paramsall,1,Wsq_diag2);
G_W2=GMM_noL(x0,pars,momentest,Wall2,momentall,paramsall,1,Wsq_all2);

W_=Wall2; % WEIGHTED!
Wsq=Wsq_all2;

global Wadd GMINadd
Wadd=Wall;
GMINadd=G_W;

%%
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
