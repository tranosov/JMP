%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation -fminsearch without L recomputing.\n');


global filename1 filename2
filename1 = "./estimation/progress_mess_econstat5.txt";
filename2 = "./estimation/progressmoment_mess_econstat5.txt";
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);

io = fopen(filename2,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);


%%
forparams =readtable('./input/PREP_econstat2.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
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

W__=diag(diag(table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames))));
Wdiag=(DOWN.*W__)\eye(size(momentest.Properties.RowNames,1));
Wsq_diag=chol(Wdiag);


%Wdiag2=Wdiag.*diag(select.('blowW'));
%Wsq_diag2=Wsq_diag.*diag(select.('blowW'));

DOWN=10^3;
Wall=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames));
Wall=(DOWN.*Wall)\eye(size(momentest.Properties.RowNames,1));
Wsq_all=chol(Wall);
scale=momentest;
scale.('blowW')=zeros(size(momentest.Properties.RowNames,1),1);
scale('scommiles','blowW')={1};
%select('swcommiles_difw','blowW')={10^3};
scale('shcommiles_dif','blowW')={1};
scale('hdj','blowW')={0.1};
scale('hdo','blowW')={0.1};
scale('sdj_dif','blowW')={0.5};
scale('wlfp_dif','blowW')={0}; %{10^3};
scale('wagegap_hw_withn','blowW')={0}; %{10^3};
scale('betahrs_w','blowW')={0.1}; %{10^2};
scale('betalwg_w','blowW')={0.1};
scale('p_gradient_simple','blowW')={0}; %{10^3};
Wall2=Wall+diag(scale.('blowW').*diag(Wall));% NOTICE - THERE WAS AN ERROR HERE. NOT ELEMENT BY ELEMNT, HAS TO BE MATRIX MULTIPLICATION. FUCK
Wsq_all2=chol(Wall2);
%Wsq_all.*diag(select.('blowW'));
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


%G_Wdiag=GMM_noL(x0,pars,momentest,Wdiag,momentall,paramsall,1,Wsq_diag);
G_W=GMM_noL(x0,pars,momentest,Wall,momentall,paramsall,1,Wsq_all);

%G_Wdiag2=GMM_noL(x0,pars,momentest,Wdiag2,momentall,paramsall,1,Wsq_diag2);
G_W2=GMM_noL(x0,pars,momentest,Wall2,momentall,paramsall,1,Wsq_all2);

W_=Wall2;
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

% options
options = optimset('Display','iter');
    ObjectiveFunction=@(x)GMM_noL(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
    rng(353);
    options.TolFun=10^(-3);
    options.TolX=50; % ignore this
    options.MaxFunctionEvaluations=3000;
    
    %Unlike other solvers, fminsearch stops when it satisfies both TolFun and TolX.
[x,fval,exitFlag,output] = fminsearch(ObjectiveFunction,x0,options);
output

% try fminunc?

io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Routine terminated. \n");
fprintf(io,"%s",output.message);
fprintf(io," \n");
fclose(io);

% 
% 1191 iterations and still did not fucking converge. still is improving.
% that is wild.
