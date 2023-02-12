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
filename1 = "./estimation/progress_sobol_econstat5.txt";
filename2 = "./estimation/progressmoment_sobol_econstat5.txt";
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
DOWN=10^3;
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Scaling down weighting matrix by a factor %16.8f\n",DOWN);
fclose(io);

Wall=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames));
Wall=(DOWN.*Wall)\eye(size(momentest.Properties.RowNames,1));
% ADD MM CLEARING CONDITION
    momentest('clmm',:)={0}; 
    Wall(end+1,end+1)=DOWN*10;

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

%G_Wdiag2=GMM_noL(x0,pars,momentest,Wdiag2,momentall,paramsall,1,Wsq_diag2);
G_W2=GMM_noL_MM(x0,pars,momentest,Wall2,momentall,paramsall,1,Wsq_all2);

W_=Wall2;
Wsq=Wsq_all2;

global Wadd GMINadd
Wadd=W_;
GMINadd=G_W2;


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
    ObjectiveFunction=@(x)GMM_noL_MM(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
    rng(354);
    options.TolFun=10^(-1); % lenient!
    options.TolX=10^(0); % hopefully  - should not matetr
    options.MaxFunctionEvaluations=500; % pretty low!

    %Unlike other solvers, fminsearch stops when it satisfies both TolFun and TolX.

ObjectiveFunction(x0)   
[x,fval,exitFlag,output] = fminsearch(ObjectiveFunction,x0,options);
output


io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Routine terminated. \n");

fprintf(io,"%s",output.message);
fprintf(io," \n");
fclose(io);




NP=size(table2array(pars),1);
sob = sobolset(NP,'Skip',1e3,'Leap',1e2);
sob= scramble(sob,'MatousekAffineOwen');
SOBN=20;
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
    fprintf(io,"%s",output.message);
    fprintf(io," \n");
    fclose(io);

    
end




% try fminunc?

% 1191 iterations and still did not fucking converge. still is improving.
% that is wild.
