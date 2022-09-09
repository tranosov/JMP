%{
Couples and singles location
Estimate GMM

%}

clear; clc;
clear global;
%cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'
rng(357)
fprintf('Running GMM estimation - F=0, fminsearch without L recomputing.\n');
global FLIN
FLIN=0;
% honestly - this is a lot I think about crrat?

% This is a bit too good for my liking - figure out why can husbands
% commute this much more?
% it has terrible participation gap, but many have, so ignore.

% mechanisms going against my point: since my wage benefit to commuting is
% proportional - it matters more for HUSBAND (it is from a bigger base)
% it seems to be such a big deal, that income of husbands is BIGGER when
% they commute!
% go back to baseline wc.

% husband hours fall by about 3 percent with a commute. if the wage
% increases by more, their income goes up with h commute.
% plus husband hours seem way less affected by commuting that wife hours!

% also - no idea why - but wife housework changes less? wife commutes-
% higher ux. does about 1/4 to 1/5 of the preference for h commute.

% this discrepency with hours is ALSO due to the pid element  - which also
% in itself makes w commuting more costly. 

% commuting single X married men - singles commute more in city, husbands
% more in suburbs.

% surprisingly, I think I also get quite a lot of commmuting couple X
% single difference from collocation problems (within each location -
% couples>singles in distance to their offer)


%AGAIN - THIS IS GOING A BIT TOO WELL???

global filename1 filename2
filename1 = "./estimation/progress_f0.txt";
filename2 = "./estimation/progressmoment_f0.txt";
io = fopen(filename1,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);

io = fopen(filename2,'a');
fprintf(io," \n");
fprintf(io,"Rerun routine. \n");
fclose(io);

%%
forparams =readtable('./input/PREP_f0.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

fprintf('Inputs loaded.\n');

% GET RID OF PHID
forparams('PHID_','toestimateR')={0};
forparams('PHID_','value')={-1000};

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

select=momentest;
select.('blowW')=ones(size(W__,1),1);
select('L','blowW')={10^3};
select('scommiles','blowW')={10^3}; %SHOULD BE A BIG CONSTRAINT ON THIS EXCERCISE!
select('swcommiles_difw','blowW')={10^2};
select('shcommiles_dif','blowW')={10^2};
select('wagegap_hw_withn','blowW')={10^2};
select('betahrs_w','blowW')={1}; % dont focus on
select('betalwg_w','blowW')={1};
select('wlfp_dif','blowW')={10^2};
Wdiag2=Wdiag.*diag(select.('blowW'));
Wsq_diag2=Wsq_diag.*diag(select.('blowW'));

DOWN=10^3;
Wall=table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames));
Wall=(DOWN.*Wall)\eye(size(momentest.Properties.RowNames,1));
Wsq_all=chol(Wall);
Wall2=Wall.*diag(select.('blowW'));
Wsq_all2=Wsq_all.*diag(select.('blowW'));
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


G_Wdiag=GMM_noL(x0,pars,momentest,Wdiag,momentall,paramsall,1,Wsq_diag);
G_W=GMM_noL(x0,pars,momentest,Wall,momentall,paramsall,1,Wsq_all);

G_Wdiag2=GMM_noL(x0,pars,momentest,Wdiag2,momentall,paramsall,1,Wsq_diag2);
G_W2=GMM_noL(x0,pars,momentest,Wall2,momentall,paramsall,1,Wsq_all2);

W_=Wall2;
Wsq=Wsq_all2;

global Wadd
Wadd=Wall;

global GMIN ITER
GMIN=99999;
ITER=0;

%%
io = fopen(filename1,'a');
fprintf(" \n");
fprintf(io," FLAG: MINIMUM \n");
fclose(io);

% options
options = optimset('Display','iter');
    ObjectiveFunction=@(x)GMM_noL(x,pars,momentest,W_,momentall,paramsall,1,Wsq); % pars has the list!
    rng(354);
    options.TolFun=10^(-1); % lenient!
    options.TolX=10^(0); % ignore this
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



%%
NP=size(table2array(pars),1);
sob = sobolset(NP,'Skip',1e3,'Leap',1e2);
sob= scramble(sob,'MatousekAffineOwen');
SOBN=10;
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
