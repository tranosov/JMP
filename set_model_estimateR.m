
function [EQS,PARREST]=set_model_estimateR(params,reinitialize)


% untransform - separate function before!


% UNCHANGED
wa=params{'wa',:};
dmin=params{'dmin',:};
b=params{'b',:};
lsbar=params{'lsbar',:};
wb=params{'wb',:};
sstaysingleh=params{'sstaysingleh',:};
sstaysinglew=params{'sstaysinglew',:};
MtoF=params{'MtoF',:};

% ESTIMATED
AS1=params{'AS1',:};
A1=params{'A1',:};
AC_sub=params{'AC_sub',:};
d=params{'d',:};
Jcenter=params{'Jcenter',:};
jobdif_=params{'jobdif_',:};
Jseg=params{'Jseg',:};
sw=params{'sw',:};
sel=params{'sel',:};
gamd_=params{'gamd_',:};
PHID_=params{'PHID_',:};
line=params{'line',:};
muw=params{'muw',:};
crra_=params{'crra_',:};
Xi=params{'Xi',:};
pid_=params{'pid_',:};
wgap_raw=params{'wgap_raw',:};
wc=params{'wc',:};
AC2=params{'AC2',:};
mm_=params{'mm',:};
sigmam=params{'sigmam',:};

LA0=params{'LA0',:};
deltaw_=params{'deltaw_',:};
plocal_=params{'plocal_',:};

% ALSO ESTIMATED
eta_=params{'eta_',:};
crrat_=params{'crrat_',:};
ce_=params{'ce_',:};


% PREPPED BEFORE
NLY=params{'NLY',:};
pish_=params{'pish_',:};
piel_=params{'piel_',:};
crrax_=params{'crrax_',:};
piw_=params{'piw_',:};
pi_=params{'pi_',:}; 

%ces_=params{'ces_',:};
%crrah_=params{'crrah_',:};


THETA=params{'THETA',:}; % but gets adjusted later!
THETAHW=params{'THETAHW',:};


%AC3=0;
%dmin=0;
dminS=0;
a=1;
ah=1; 
aw=1;
bh=1;
bw=aw; 
timeh=0; 
timew=timeh; 
Jschool=1;
Wcenter=1;
Hcenter=1;
wlinear=2;
cneces=0;
mI_=0; 
mA_=0;
kappa_=0;
PHI_=0;
PHIW_=0; 
rhod_=0; 
pitheta_=0;
qel_=0;
d21=1;
d31=1;
Sup1=0;
Sup2=0;
Sup3=0;
SupS=0;
addS=0;
out=1;
daddS=1;
aaddS=1;
pophous=1; 
%wc=0;
PID=0; 
%lsbar=1;

ces_=ce_; % fixed
crrah_=crra_;


ceh_=ce_;
cew_=ce_;
leh=1;
lew=1;
etaC_=eta_; % fixed



Wschool=1/jobdif_; % close to 1 - evenly distributed.
Hschool=jobdif_; 
muprob=1;
muwL=muw;
swL=sw;
pih_=1-piw_;
pisw_=pish_;

typeic_=3.5;

toinputs=reinitialize;
[EQS,PARREST]=...
    model_new(toinputs,AS1,A1,AC_sub, AC2,dmin,dminS,d, a,ah,aw,b,bh,bw,timeh,timew,LA0,Jcenter, Jschool,...
    Wcenter, Wschool,Hcenter, Hschool,Jseg,sw,sel,wa,wb,wlinear, cneces,...
    ces_,ceh_,cew_,leh,lew,eta_,etaC_, muprob,plocal_,sstaysingleh,sstaysinglew,mA_,mI_,kappa_,...
    gamd_,PHID_,deltaw_,PHI_,rhod_,PHIW_,line,muw,muwL,swL,crra_,crrat_,crrah_,crrax_,pih_,piw_,pi_,pish_,pisw_,piel_,pitheta_,qel_,...
    d21,d31,Sup1,Sup2,Sup3,SupS,addS,out,daddS,aaddS,pophous, Xi,wc, PID,typeic_,pid_,wgap_raw,mm_,lsbar,NLY,THETA,THETAHW,sigmam,MtoF);


params('crrah_','value')={'crrah_'};
params('ces_','value')={'ces_'};
PARREST.('params')=params;


end
