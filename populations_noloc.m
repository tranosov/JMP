
function [DC, DC1, DC2,HCouple,Pw,hC,VC_,cexp1,cexp2,...
    OUTC, DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS,EQS,PARREST,hS] = populations_noloc(p,LA,EQS,PARREST)

global VERBOSE doepses doublecouple
if isempty(doepses)
    doepses=0;
end

VC_=99; %what is this?
D=PARREST.('D');
JLs=PARREST.('JLs');
NC=PARREST.('NC');
mu=PARREST.('mu');
NS=PARREST.('NS');
NSh=PARREST.('NSh');
NSw=PARREST.('NSw');

%DC=PARREST.('DC');
DC1=PARREST.('DC1');
%DSh=PARREST.('DSh');
%DSw=PARREST.('DSw');
DS=PARREST.('DS');

wfh=PARREST.('wfh');
HS=PARREST.('HS');
I=size(HS,2);
T=size(JLs,2);
if wfh>0
    T=T*2;
end
W=3;

[~, ~, ~,~,Pw,hC,~,cexp1,cexp2,OUTC, ~,~, ~, ~, VS,Pws,cexps,OUTS,EQS,PARREST,hS] = ...
        populations(p,LA,EQS,PARREST,0); 
    

PS=sum(DS,4)./repmat(sum(sum(DS,4),3),1,1,3);
PS(isnan(PS))=0;
DS=Pws.*repmat(NS,1,1,I,W).*repmat(PS,1,1,1,W);
DSh=Pws.*repmat(NSh,1,1,I,W).*repmat(PS,1,1,1,W);
DSw=Pws.*repmat(NSw,1,1,I,W).*repmat(PS,1,1,1,W);
HSingleH=hS.*DSh; % no point recomputingx share single here
HSingleW=hS.*DSw;
HSingle=HSingleH+HSingleW;

PC=sum(sum(DC1,6),7)./repmat(sum(sum(sum(DC1,6),7),5),1,1,1,1,3);
PC(isnan(PC))=0;
DC1=Pw.*repmat(PC,1,1,1,1,1,W,W).*repmat(NC,1,1,1,1,I,W,W)./repmat(2,I,I,T,T,I,W,W);
pjt=NC./repmat(sum(sum(NC)),I,I);
pjt(isnan(pjt))=0;

if doublecouple
    DC2=Pw.*((mu^2)*repmat(pjt,1,1,1,1,I,W,W).*repmat(sum(sum(sum(sum(DC1,7),6))),I,I,1,1,1,W,W)+...
    mu*(1-mu)*repmat(repmat(sum(pjt,2),1,I,1,1),1,1,1,1,I,W,W).*repmat(sum(sum(sum(DC1,6),7)),I,1,1,1,1,W,W)+...
    mu*(1-mu)*repmat(repmat(sum(pjt,1),I,1,1,1),1,1,1,1,I,W,W).*repmat(sum(sum(sum(DC1,6),7),2),1,I,1,1,1,W,W)+...
    ((1-mu)^2)*repmat(sum(sum(DC1,6),7),1,1,1,1,1,W,W)); 
else
    DC2=DC1;
end

DC=DC1+DC2;

HCouple=DC.*hC;

% is it ok that there is not epses here?
end