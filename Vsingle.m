%{
Couples and singles location

Function computing the value for singles.

t - type
j - job offer location
i - residence location
p - price of housing

Maybe later add wage, now treat as a global variable!
%}

function [VS,works,Pn,Pw0,Pw,conexp]=Vsingle(t,j,i, PARREST,OUTS)

plocal=PARREST.('plocal');

%[mIL_,mAL_,mI_,mA_]=matchdist(i,j,t,mA,mI,mAL,mIL, typeic,D);
%mA_-mAL_
[VS_,works_,Pn_,Pw0_,Pw_,conexp_]=Vsingle0(t,j,i,PARREST,OUTS);
[VS__,works__,Pn__,Pw0__,Pw__,conexp__]=Vsingle1(t,j,i,PARREST,OUTS);

output=plocal.*[VS__,works__,Pn__,Pw0__,Pw__,conexp__]+...
    (1-plocal).*[VS_,works_,Pn_,Pw0_,Pw_,conexp_];
VS=output(1);
works=output(2);
Pn=output(3);
Pw0=output(4);
Pw=output(5);
conexp=output(6);
%Pw__ 
%I do not think mean of dist matters. increasing sigma increases this, it gets acceptence closer to a half

end
