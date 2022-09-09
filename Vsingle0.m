%{
Couples and singles location

Function computing the value for singles.

t - type
j - job offer location
i - residence location
p - price of housing

Maybe later add wage, now treat as a global variable!
%}

function [VS,works,Pn,Pw0,Pw,conexp]=Vsingle0(t,j,i,PARREST,OUTS)

vs=OUTS.('vs') ;
%typeic=PARREST.('typeic');
%D=PARREST.('D');
%mm=PARREST.('mm');
%JLs=PARREST.('JLs');
%betah=PARREST.('betah');

mIL=PARREST.('mIL_');
mAL=PARREST.('mAL_');
mI=PARREST.('mI_');
mA=PARREST.('mA_');

mIL_=mIL(t,j,i);
mAL_=mAL(t,j,i);
mI_=mI(t,j,i);
mA_=mA(t,j,i);

%[mIL_,mAL_,mI_,mA_]=matchdist(i,j,t,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);

vw=vs(t,j,i); %us(1,d,a,cs(Yss(t,j,i),p),hdS(Yss(t,j,i),p),lss(t,j,i),ics(t,j,i)); 
%vw0=vs(t,i,i,p); %us(1,d0,a,cs(Yss(t,i,i),p),hdS(Yss(t,i,i),p),lss(t,i,i),ics(t,i,i)); 
vn=-1.0000e+10; %us(0,0,a,0,0,0,0);
Pw=1;
Pn=0;
Pw0=0;
works=Pw; 
conexp=(mI_+mA_)/2;
VS=(1-Pw)*vn+Pw*vw+conexp;

end

% notice - the sequential nature here basically means that the local option 
% is treated as the mean - its dispersion does not matter!
