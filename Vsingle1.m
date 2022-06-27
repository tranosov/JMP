%{
Couples and singles location

Function computing the value for singles.

t - type
j - job offer location
i - residence location
p - price of housing

Maybe later add wage, now treat as a global variable!
%}

function [VS,works,Pn,Pw0,Pw,conexp]=Vsingle1(t,j,i,PARREST,OUTS)
global show
vs=OUTS.('vs') ;
typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
betah=PARREST.('betah');
mA=PARREST.('mA');
mI=PARREST.('mI');
mAL=PARREST.('mAL');
mIL=PARREST.('mIL');

[mIL_,mAL_,mI_,mA_]=matchdist(i,j,t,mA,mI,mAL,mIL,typeic,D,mm,JLs, betah);

vw=vs(t,j,i); %us(1,d,a,cs(Yss(t,j,i),p),hdS(Yss(t,j,i),p),lss(t,j,i),ics(t,j,i)); 
vw0=vs(t,i,i); %us(1,d0,a,cs(Yss(t,i,i),p),hdS(Yss(t,i,i),p),lss(t,i,i),ics(t,i,i)); 
vn=-1.0000e+10; %us(0,0,a,0,0,0,0);
[vbar,vbari]=max([vn,vw]);
conexp_=(mA_+mI_)/2; %between vw and not working - always pick vw
vbar=vbar+conexp_; % todo: change this (and the couples integrals) to deciding simultaneously??
Pw0= (mAL_-max(mIL_,min(mAL_,(vbar-vw0))))/(mAL_-mIL_);
if vbari==1
    Pw=0;
    Pn=1-Pw0; % this will not come up, so ignore.
else
    Pw=1-Pw0;
    Pn=0; 

works=Pw+Pw0; % ultimately I think this will NOT  be useful
conexp=Pw0*(mA_+max(mIL_,min(mAL_,(vbar-vw0))))/2;
VS=Pw*vbar+Pw0*vw0+conexp;
conexp=Pw*conexp_+ conexp;
if show==1
vw0-vw    
Pw=Pw %- for this to be enough need to be around 37-42 for a third of them getting an offer
end
end

% notice - the sequential nature here basically means that the local option 
% is treated as the mean - its dispersion does not matter!

%order local-global does not matter, at least for Pw!

%0.3 difference - without any home function!