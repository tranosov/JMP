%{
Couples and singles location

Function computing the value for couples in the second period.

t - type (does not matter here)
j - job offer location
i - residence location
p - price of housing

Maybe later add wage, now treat as a global variable!
%}



function [V,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2]=...
    Vcouplesharesmatch_2(th,tw,jh,jw,i, PARREST,OUTC)
%global plocal D  %AS  AC w1 w2 xh xw checklsb crra uw uh us ch cw cs hdS hdC Ys Yc lss_G lsah_G lsaw_G lsb_G alphaw betaw alphah betah tih tiw mu JLs Jw Jm HS Jc NC NS NSh NSw sigmaw sigmal mA mI kappa PHI mAL mIL
%global typeic mm show JLs mA mI mAL mIL OUTC
global show

typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
betah=PARREST.('betah');
mA=PARREST.('mA');
mI=PARREST.('mI');
mAL=PARREST.('mAL');
mIL=PARREST.('mIL');

plocal=PARREST.('plocal');

[mI1_,mA1_,mIG1_,mAG1_]=matchdist(i,jh,th,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);
[mI2_,mA2_,mIG2_,mAG2_]=matchdist(i,jw,tw,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);


[V_,workh_, workw_, Pnn_, Pnw0_,Pnw_, Pw0n_,Pw0w0_,Pw0w_, Pwn_, Pww0_,Pww_,conexp1_, conexp2_]=Vcouplesharesmatch0_2(th,tw,jh,jw,i,PARREST,OUTC);
[V__,workh__, workw__, Pnn__, Pnw0__,Pnw__, Pw0n__,Pw0w0__,Pw0w__, Pwn__, Pww0__,Pww__,conexp1__, conexp2__]=Vcouplesharesmatch1_2(th,tw,jh,jw,i,PARREST,OUTC);
%Pww0__-Pw0w__
V=(1-plocal).*V_+plocal*V__;
workh=(1-plocal).*workh_+plocal*workh__;
workw=(1-plocal).*workw_+plocal*workw__;
Pnn=(1-plocal).*Pnn_+plocal*Pnn__;
Pnw0=(1-plocal).*Pnw0_+plocal*Pnw0__;
Pnw=(1-plocal).*Pnw_+plocal*Pnw__;
Pw0n=(1-plocal).*Pw0n_+plocal*Pw0n__;
Pwn=(1-plocal).*Pwn_+plocal*Pwn__;
Pww=(1-plocal).*Pww_+plocal*Pww__;
Pww0=(1-plocal).*Pww0_+plocal*Pww0__;
Pw0w=(1-plocal).*Pw0w_+plocal*Pw0w__;
Pw0w0=(1-plocal).*Pw0w0_+plocal*Pw0w0__;
conexp1=(1-plocal).*conexp1_+plocal*conexp1__;
conexp2=(1-plocal).*conexp2_+plocal*conexp2__;

%{
Pww__
Pww0__
Pw0w__
Pwn__
plocal
Pwn+Pww+Pww0
Pwn__+Pww__+Pww0__
Pnw+Pww+Pw0w
Pnw__+Pww__+Pw0w__
workh
workw
%}
check=round(Pww+Pwn+Pw0w0+Pw0n+Pw0w+Pww0+Pnw+Pnw0,2);
if (check==1)*(Pww<=1)...
~=1
    check
    Pww<=1
    return
end
%{
Pww
Pw0w0
Pww0
Pw0w
min(dh,dw)

difv=vwbar-vbarw %Vbarw<Vwbar
difbar=Pwbar-Pbarw
dif=Pww0-Pw0w
Pvbar
Pwn
Pw0n
Pww
Pw0w0
Pw0w
Pww0



doh=Do(i,JLs(:,th))
dow=Do(i,JLs(:,tw))
meanmatchh= (mA1_+mI1_)/2
meanmatchw= (mA2_+mI2_)/2
meanmatchLh=(mAG1_+mIG1_)/2
meanmatchLw=(mAG2_+mIG2_)/2
hcom=Pww0+Pwn+Pww
wcom=Pw0w+Pnw+Pww
%}

%Pwn-Pnw
%Pw0n-Pnw0
%vwn-vnw
%vw0n-vnw0
if show==1
workh 
workw
end
end
