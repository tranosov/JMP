%{
Couples and singles location

Function computing the value for couples in the second period.

t - type (does not matter here)
j - job offer location
i - residence location
p - price of housing

Maybe later add wage, now treat as a global variable!

% LOCAL FIRST
%}



function [V,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2]=...
    Vcouplesharesmatch1_2(th,tw,jh,jw,i, PARREST,OUTC)
%global D  %AS  AC w1 w2 xh xw checklsb crra uw uh us ch cw cs hdS hdC Ys Yc lss_G lsah_G lsaw_G lsb_G alphaw betaw alphah betah tih tiw mu  JLs Jw Jm HS Jc NC NS NSh NSw sigmaw sigmal mA mI kappa PHI mAL mIL
%global OUTC typeic show EQS mm mA mI mAL mIL %lsah lsaw lsb Ycc ich icw
global show

mIL=PARREST.('mIL_');
mAL=PARREST.('mAL_');
mI=PARREST.('mI_');
mA=PARREST.('mA_');

mI1_=mIL(th,jh,i);
mA1_=mAL(th,jh,i);
mIG1_=mI(th,jh,i);
mAG1_=mA(th,jh,i);
mI2_=mIL(tw,jw,i);
mA2_=mAL(tw,jw,i);
mIG2_=mI(tw,jw,i);
mAG2_=mA(tw,jw,i);


lambda=OUTC.('lambda'); %this way I am sure it was used in vc!
vc=OUTC.('vc') ;
 
% choice: not working is 0 commute!

%todo: add types to lsb lsah lsaw
% ORDER: HUSBAND-WIFE
vw0w0=vc(th,tw,i,i,i,3);
vw0n=vc(th,tw,i,i,i,1);
vnw0=vc(th,tw,i,i,i,2);
vww=vc(th,tw,jh,jw,i,3);
vwn=vc(th,tw,jh,jw,i,1);
vnw=vc(th,tw,jh,jw,i,2);
vnn=-1.0000e+10;
vww0=vc(th,tw,jh,i,i,3);
vw0w=vc(th,tw,i,jw,i,3);

if vnn>=max([vnw0,vw0n,vw0w0])
    vww 
    vwn 
    vnw 
    vnn
    return
end


[V,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2]=...
        locfirst(i,jh,jw,th,tw,lambda,mI1_,mA1_,mIG1_,mAG1_,mI2_,mA2_,mIG2_,mAG2_,vw0w0,vw0n,vnw0,vww,vwn,vnw,vnn,vww0,vw0w);
%{
workh
workw
Pww0-Pw0w
Pww+Pww0+Pwn
Pww+Pw0w+Pnw
Pww
%Pww0+Pwn
%Pw0w+Pnw
%Pw0w0
[V,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2]=globfirst(i,jh,jw,th,tw,lambda,mA,mI,mAL,mIL,typeic,vw0w0,vw0n,vnw0,vww,vwn,vnw,vnn,vww0,vw0w);


%}
%Pww0-Pw0w
%workh-workw
%Pwn-Pnw
%Pw0n-Pnw0
%vwn-vnw
%vw0n-vnw0
if show==1
    %vww0
    %vw0w0
    dw_hus=-vw0w0+vww0
    dw_wif=-vw0w0+vw0w
%workh
%workw
%Pww0-Pw0w
haccepts=Pww+Pww0+Pwn
waccepts=Pww+Pw0w+Pnw
%Pww
Pw0w0

vww-vwn
vww-vnw
end
%Pww0+Pwn
%Pw0w+Pnw
%Pw0w0
%}
end


% btw - order does not seem to matter if there is not participation, but it
% does if there is.
% loc first - supports the notion that selection -> more men commute. the
% other one does not??

% dh,dw incrase - slightly icnreases the difference between vwn and vnw
% notice - if one at home - commuting of other one essentially does NOT
% matter with min(dh,dw), no matter what delta, because it is still 0!


% issue underlining now - without additional help couples commute LESS -
% why? 0.18 is the welfare loss of h commuting. for singles it is bigger!
% (0.3!)

% I am not at all sure why, but - increasing dispersion makes a lot of
% singles accept more offers, but for couples increases just a tiny bit??
% and this depends a ton on the 'elasticity'. Makes the couples advantage
% potent.
% je ten increase se sigmaw nejak nonlinear?


% if I make the other option enticing - still higher dispersion makes
% singles accept more? This I positively do not get. I thought I would have
% flipped the situation!
% there is some issue I am NOT getting. The thing is bridging the gap with
% match shock - couples basically get just a half of it, because almost
% always just 1 commutes.
% with my function - I need to be able to make it so that for singles ( am
% bellow 50% and for at least one commuting I am above!

% 14-20-22 - 0.6
% 10.7-18-20 - 0.55.

