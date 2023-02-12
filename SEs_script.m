clear; clc;
clear global;
cd 'C:\Users\Admin\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'

rng(357)
global show %mm typeic OUTC OUTS  sigmam params betah %mA mI kappa PHI mAL mIL ceh cew eta etaC XI plocal pih piw PHID Fdist FdistS D  AS  AC w1 w2 xh xw checklsb crra uu uw uh us ch cw cs hdS hdC Ys Yc lss_G lsah_G lsaw_G lsb_G alphaw betaw alphah betah tih tiw mu lambda JLs Jw Jm HS Jc NC NS NSh NSw sigmaw sigmal 
show=0;
global VERBOSE doepses
VERBOSE=1;
doepses=1;

toplot="C:\Users\Admin\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab\output\";
addon="d"; %econstat2"; 
fileall =toplot+"SEs"+'_'+addon+'.xlsx';

%%
forparams =readtable('./input/PREP_' +addon + '.xlsx','Sheet','PARS','ReadVariableNames', true,'ReadRowNames',true);
formom =readtable('./input/PREP.xlsx','Sheet','MOMS','ReadVariableNames', true,'ReadRowNames',true);
forw =readtable('./input/PREP.xlsx','Sheet','W','ReadVariableNames', true,'ReadRowNames',true);

global FLIN
FLIN=0;
params=forparams(:,'value');
paramsest=forparams(forparams.('toestimateR')==1,'value');
pars=paramsest;
pars_=table2array(paramsest);
params(pars.Properties.RowNames,:)= array2table(pars_, 'RowNames',pars.Properties.RowNames);
params=untransform(params);

formom('do_sd','toestimate')={1} % ADD BACK TO SEE
momentest=formom(formom.('toestimate')==1,'value');
parsc=calibrate1(params,momentest); %calibrate within
params(parsc.Properties.RowNames,:)=parsc;
[EQS,PARREST] = set_model_estimateR(params,1);
params=PARREST.('params');

p0=solvemodel(params,EQS,PARREST,[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ],params{'LA0',:},0,1,0);

[DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,...
    OUTC, DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS,EQS,PARREST] = ...
        populations(p0,params{'LA0',:},EQS,PARREST,0.1);  %recalibrates thetas! % and before add E(eps)!
params=PARREST.('params');

params('p0_1','value')={p0(1)};
params('p0_2','value')={p0(2)};
params('p0_3','value')={p0(3)};


COV=(table2array(forw( momentest.Properties.RowNames, momentest.Properties.RowNames)));
W=COV\eye(size(COV,1));
W_diag=diag(diag(COV))\eye(size(COV,1));
%array2table(diag(W)-diag(W_diag),'RowNames',momentest.Properties.RowNames)
%Why are diagonals on W much bigger than W diag?
momentall=formom(:,'value');


pars0=params(forparams.('forses')==1,'value');

%%
pars=pars0;
selectp=pars;
selectp.('forses')=ones(size(table2array(pars0),1),1);
selectp('LA0','forses')={0};
pars=pars(selectp.('forses')==1,:);
pars_=table2array(pars); % 27!

n=size(table2array(momentest),1);
noL=0;
n_=n;

jac0sheet='jac0';
if isfile(fileall)
    [~,sheets] = xlsfinfo(fileall);
    if any(strcmp(sheets,jac0sheet))
        jac0_=table2array(readtable(toplot+"SEs"+'_'+addon+'.xlsx','Sheet',jac0sheet,'ReadVariableNames', false,'ReadRowNames',false));
        jac0_=jac0_(1:n_,:);
    else 
         jac0_=zeros(n_,size(pars_,1));
    end
else 
    jac0_=zeros(n_,size(pars_,1));
end
iii=1;
%[jac0]=SEs(pars_,pars,momentest,W,momentall,params,fileall,jac0sheet,n_,iii,jac0_,0,0)

[jac0]=SEs(pars_,pars,momentest,W,momentall,params,fileall,jac0sheet,n_,iii,jac0_,0,0)
% 0 recalibrate - do not recalibrate any parameters, and not even thetas
% 0 withmm - always assume the moment with lambda is lambda (do not
% computethe numeric derivative!)

%{
pars=pars0;
selectp=pars;
selectp.('forses')=ones(size(table2array(pars0),1),1);
selectp('THETAHW','forses')={0};
pars=pars(selectp.('forses')==1,:);
pars_=table2array(pars);

n=size(table2array(momentest),1);
noL=1;
n_=n;
if noL
    n_=n-1;% take 1 moment out
end
sh='jac0_lambda';
if isfile(fileall)
    [~,sheets] = xlsfinfo(fileall);
    if any(strcmp(sheets,sh))
        jac0_=table2array(readtable(toplot+"SEs"+'_'+addon+'.xlsx','Sheet',sh,'ReadVariableNames', false,'ReadRowNames',false));
        jac0_=jac0_(1:n_,:);
    else 
         jac0_=zeros(n_,size(pars_,1));
    end
else 
    jac0_=zeros(n_,size(pars_,1));
end
iii=1;
[jac0_lambda]=SEs(pars_,pars,momentest,W,momentall,params,fileall,jac0sheet,n_,iii,jac0_,0.5,1-noL)
%}

%%
jac0sheet='jac0_8_2'; % better than jac0
addon='d';
jac0=table2array(readtable(toplot+"SEs"+'_'+addon+'.xlsx','Sheet',jac0sheet,'ReadVariableNames', false,'ReadRowNames',false));
select=momentest;
select.('forses')=ones(size(jac0,1),1);
select('NLY_','forses')={0};
select('shouseexp','forses')={0};
select('L','forses')={1}; % why not? I would say keep it?

pars=pars0;
pars_=table2array(pars);
selectp=pars;
selectp.('forses')=ones(size(jac0,2),1);
%selectp('ce_','forses')={1}; % is just one in the new version
selectp('eta_','forses')={0};
%selectp('THETAHW','forses')={0};
jac0_v2=jac0(logical(table2array(select(:,'forses'))),logical(table2array(selectp(:,'forses'))));
G0=jac0_v2';

COV=(table2array(forw(momentest(select.('forses')==1,:).Properties.RowNames, momentest(select.('forses')==1,:).Properties.RowNames)));
W=COV\eye(size(COV,1)); % NO WEIGHTING HERE
W_diag=diag(diag(COV))\eye(size(COV,1));

VAR_=G0*W*G0';
VAR=VAR_\eye(size(VAR_,1));
out0=array2table([sqrt( diag(VAR)), pars_(logical(table2array(selectp(:,'forses'))))./sqrt( diag(VAR))],'RowNames',selectp(selectp.('forses')==1,'forses').Properties.RowNames,'VariableNames',{'se','t'})


writetable(out0,fileall,'Sheet','SEs1','WriteRowNames',true) % MAIN - in the paper
%%
jac0sheet='jac0_8_2'; % better than jac0

jac0=table2array(readtable(toplot+"SEs"+'_'+addon+'.xlsx','Sheet',jac0sheet,'ReadVariableNames', false,'ReadRowNames',false));
select=momentest;
select.('forses')=ones(size(jac0,1),1);
select('NLY_','forses')={0};
select('shouseexp','forses')={0};
select('L','forses')={1};

pars=pars0;
pars_=table2array(pars);
selectp=pars;
selectp.('forses')=ones(size(jac0,2),1);
selectp('ce_','forses')={1};
selectp('eta_','forses')={0};
selectp('THETAHW','forses')={1};
jac0_v2=jac0(logical(table2array(select(:,'forses'))),logical(table2array(selectp(:,'forses'))));
G0=jac0_v2';

COV=(table2array(forw(momentest(select.('forses')==1,:).Properties.RowNames, momentest(select.('forses')==1,:).Properties.RowNames)));
W=COV\eye(size(COV,1));
W_diag=diag(diag(COV))\eye(size(COV,1));

VAR_=G0*W*G0';
VAR=VAR_\eye(size(VAR_,1));
out0=array2table([sqrt( diag(VAR)), pars_(logical(table2array(selectp(:,'forses'))))./sqrt( diag(VAR))],'RowNames',selectp(selectp.('forses')==1,'forses').Properties.RowNames,'VariableNames',{'se','t'})


writetable(out0,fileall,'Sheet','SEs1','WriteRowNames',true) % MAIN- but no

pars=pars0;
pars_=table2array(pars);
selectp=pars;
selectp.('forses')=ones(size(jac0,2),1);
selectp('ce_','forses')={1};
selectp('eta_','forses')={0};
selectp('THETAHW','forses')={0};
jac0_v2=jac0(logical(table2array(select(:,'forses'))),logical(table2array(selectp(:,'forses'))));
G0=jac0_v2';

COV=(table2array(forw(momentest(select.('forses')==1,:).Properties.RowNames, momentest(select.('forses')==1,:).Properties.RowNames)));
W=COV\eye(size(COV,1));
W_diag=diag(diag(COV))\eye(size(COV,1));

VAR_=G0*W*G0';
VAR=VAR_\eye(size(VAR_,1));
out0=array2table([sqrt( diag(VAR)), pars_(logical(table2array(selectp(:,'forses'))))./sqrt( diag(VAR))],'RowNames',selectp(selectp.('forses')==1,'forses').Properties.RowNames,'VariableNames',{'se','t'})


writetable(out0,fileall,'Sheet','SEs1_noTHETAHW','WriteRowNames',true) % MAIN - in the paper, not anymore


% getting rid of THETAHW helps a lot. I think that was my issure to be
% honest...
% this is true EVEN WHEN the L moment is kicked out. 

% GETTING RID OF pid DOES NOTHING for me...
% also true for deltaw_! in fact this seems to maybe make things worse??
% only THETA would also help I guess. but then I think I would add sem and
% maybe would be in the same place? (PHID is not great but is usable)

% wouldn't it be than more honest to go to a version of without resolving
% for lambda!? and treat lambda as calibrated outside? Sure, maybe, but
% then I would not have market clearing in the marriage market when
% computing these derivatives? that also does not make sense to me.

DOWN=10^3;
select=momentest;
select.('forses')=ones(size(jac0,1),1);
select('NLY_','forses')={0};
select('shouseexp','forses')={0};
select('L','forses')={1};
Wall=(table2array(forw(momentest(select.('forses')==1,:).Properties.RowNames, momentest(select.('forses')==1,:).Properties.RowNames)));
Wall=(DOWN.*Wall)\eye(size(momentest(select.('forses')==1,:).Properties.RowNames,1));

scale=momentest(select.('forses')==1,:);
scale.('blowW')=zeros(size(momentest(select.('forses')==1,:).Properties.RowNames,1),1);
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
Ww=Wall+diag(scale.('blowW').*diag(Wall));

VARin=(G0*Ww*G0')\eye(size(VAR,1));
VAR=VARin*(G0*Ww*COV*Ww'*G0')*VARin;
out=array2table([sqrt( diag(VAR)), pars_(logical(table2array(selectp(:,'forses'))))./sqrt( diag(VAR))],'RowNames',selectp(selectp.('forses')==1,'forses').Properties.RowNames,'VariableNames',{'se','t'})

% This is PRETTY SENSITIVE, BUT  a little bit on the commuting moments I
% can work with. It is the betas that makes this worse, in a way seems
% unavoidable. Same with hdj or hdo!

writetable(out,fileall,'Sheet','SEs1_noTHETAHW_weighted','WriteRowNames',true)

%%
VARin=(G0*W_diag*G0')\eye(size(VAR,1));
VAR=VARin*(G0*W_diag*COV*W_diag'*G0')*VARin;
out=array2table([sqrt( diag(VAR)), pars_(logical(table2array(selectp(:,'forses'))))./sqrt( diag(VAR))],'RowNames',selectp(selectp.('forses')==1,'forses').Properties.RowNames,'VariableNames',{'se','t'})
writetable(out,fileall,'Sheet','SEs1_wdiag','WriteRowNames',true)


%%
pars=pars0;
selectp=pars;
selectp.('forses')=ones(size(table2array(pars0),1),1);
selectp('LA0','forses')={0};
pars=pars(selectp.('forses')==1,:);
pars_=table2array(pars);

jac0sheet='jac0';
jac0=table2array(readtable(toplot+"SEs"+'_'+addon+'.xlsx','Sheet',jac0sheet,'ReadVariableNames', false,'ReadRowNames',false));
G0=jac0' ;
VAR=G0*W*G0';
VAR=VAR\eye(size(VAR,1));
out0=array2table([sqrt( diag(VAR)),pars_./sqrt( diag(VAR))],'RowNames',pars.Properties.RowNames,'VariableNames',{'se','t'})
writetable(out0,fileall,'Sheet','SEs0','WriteRowNames',true) 

VARin=(G0*W_diag*G0')\eye(size(VAR,1));
VAR=VARin*(G0*W_diag*COV*W_diag'*G0')*VARin;
out=array2table([sqrt( diag(VAR)),pars_./sqrt( diag(VAR))],'RowNames',pars.Properties.RowNames,'VariableNames',{'se','t'})
writetable(out,fileall,'Sheet','SEs0_wdiag','WriteRowNames',true)



%%
% interesting how this works - trying to compute ses with eta as well - and
% not having shouseexp - results in blowing up ses for ce and A1. I guess
% because these parameters affect moments similarly, except for shouseexp! 
% so this method is smarter than I though!

select=momentest;
select.('forses')=zeros(size(jac0,1),1);
select('NLY_','forses')={0};
select('shouseexp','forses')={1};
select('shours','forses')={1};
select('shwk','forses')={1};
select('p_gradient','forses')={0};
select('p_gradient_simple','forses')={0};

selectp=pars;
selectp.('forses')=zeros(size(jac0,2),1);
selectp('ce_','forses')={1};
selectp('eta_','forses')={0};
selectp('A1','forses')={0};
jac0_v2=jac0(logical(table2array(select(:,'forses'))),logical(table2array(selectp(:,'forses'))));
G0=jac0_v2';

COV=(table2array(forw(momentest(select.('forses')==1,:).Properties.RowNames, momentest(select.('forses')==1,:).Properties.RowNames)));
W=COV\eye(size(COV,1));
W_diag=diag(diag(COV))\eye(size(COV,1));

VAR_=G0*W*G0';
VAR=VAR_\eye(size(VAR_,1));
out=array2table([sqrt( diag(VAR)), pars_(logical(table2array(selectp(:,'forses'))))./sqrt( diag(VAR))],'RowNames',selectp(selectp.('forses')==1,'forses').Properties.RowNames,'VariableNames',{'se','t'})

VARin=(G0*W_diag*G0')\eye(size(VAR,1));
VAR=VARin*(G0*W_diag*G0'*COV*W_diag'*G0')*VARin;
out=array2table([sqrt( diag(VAR)),pars_./sqrt( diag(VAR))],'RowNames',pars.Properties.RowNames,'VariableNames',{'se','t'})
