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
    Vcouplesharesmatch0_2(th,tw,jh,jw,i, PARREST,OUTC)
%global D  AS  AC w1 w2 xh xw checklsb crra uw uh us ch cw cs hdS hdC Ys Yc lss_G lsah_G lsaw_G lsb_G alphaw betaw alphah betah tih tiw mu  JLs Jw Jm HS Jc NC NS NSh NSw sigmaw sigmal mA mI kappa PHI mAL mIL
%global OUTC typeic EQS mm %lsah lsaw lsb Ycc ich icw

typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
betah=PARREST.('betah');
mA=PARREST.('mA');
mI=PARREST.('mI');
mAL=PARREST.('mAL');
mIL=PARREST.('mIL');


lambda=OUTC.('lambda'); 
vc=OUTC.('vc') ;

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

%[mI1_,mA1_,mIG1_,mAG1_]=matchdist(i,jh,th,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);
%[mI2_,mA2_,mIG2_,mAG2_]=matchdist(i,jw,tw,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);

 
% choice: not working is 0 commute!

%todo: add types to lsb lsah lsaw
% ORDER: HUSBAND-WIFE
vw0w0=vc(th,tw,i,i,i,3);
vw0n=vc(th,tw,i,i,i,1);
vnw0=vc(th,tw,i,i,i,2);
vww=vc(th,tw,jh,jw,i,3);
vwn=vc(th,tw,jh,i,i,1);
vnw=vc(th,tw,jh,jw,i,2);
vnn=-1.0000e+10;
vww0=vc(th,tw,jh,i,i,3);
vw0w=vc(th,tw,i,jw,i,3);
% not used!

if vnn>=max([vnw0,vw0n,vw0w0])
    vww 
    vwn 
    vnw 
    vnn
    return
end

denomG=(mAG1_-mIG1_)*(mAG2_-mIG2_);

if denomG>0
    db_x=mIG2_*(1-lambda)/lambda - (vwn-vnw)/(lambda); %line intercept
    db_y=mIG1_*(lambda)/(1-lambda) + (vwn-vnw)/((1-lambda));
    db_X=mAG2_*(1-lambda)/lambda - (vwn-vnw)/(lambda);
    db_Y=mAG1_*(lambda)/(1-lambda) + (vwn-vnw)/((1-lambda));
%not working or taking a global offer if partner also has a global offer   
    tresh1=max(min([db_X,(vnw-vww)/(lambda),mAG1_]),mIG1_);
    tresh2=max(min([db_Y,(vwn-vww)/((1-lambda)),mAG2_]),mIG2_);
    tresh1_=min(max([db_x,mIG2_*(1-lambda)/lambda+(vnw-vwn)/((lambda)),mIG1_]),mAG1_); %todo: simplify- this seems like db_X
    tresh2_=min(max([db_y,mIG1_*(lambda)/(1-lambda)+(vwn-vnw)/((1-lambda)),mIG2_]),mAG2_);
    
    MM_=max(0,min(mAG2_-mIG2_,mAG2_-(vwn-vww)/(1-lambda)))...
        *max(0,min(mAG1_-mIG1_,mAG1_-(vnw-vww)/(lambda)))/denomG; 
    Pww_=MM_;
    NN_=max(0,mAG2_-max([mIG2_,tresh2]))*...
        max(0,min(mAG1_-mIG1_,tresh1 - mIG1_ ))/denomG;
    OO_=max(0,min(mAG2_-mIG2_,tresh2 - mIG2_))*...
        max(min(mAG1_-mIG1_,mAG1_-tresh1  ),0)/denomG;
    PP_=0.5*max(0,...
        tresh2-tresh2_)*...
        max(0,tresh1-tresh1_)...
        *(MM_+NN_+OO_<1)/denomG; % does this cover all the zero cases?
    RR_=max(0,...
        min(tresh2_-mIG2_,mAG2_-mIG2_))*...
        max(0,min(tresh1-mIG1_,mAG1_-mIG1_))...
        *(MM_+NN_+OO_<1)*(mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))>=mIG2_)/denomG...
        +...
        max(0,...
        min(tresh2-mIG2_,mAG2_-mIG2_))*...
        max(0,(tresh1_-mIG1_)...
        )*(MM_+NN_+OO_<1)*(mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))<mIG2_)/denomG;
    
    Pwn_=OO_+PP_+RR_*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))>=mIG2_);
    Pnw_=NN_+PP_+RR_*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))<mIG2_);
    
    %todo: replace tresh1 when appropriate
    conexp1_=(mAG1_+tresh1)*0.5*(OO_+MM_) +...
	     (tresh1+mIG1_)*0.5*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/(1-lambda)>=mIG2_)*RR_+...
         triangle(tresh1_,tresh1,tresh2_,tresh2,0,0)*PP_;
        
    conexp2_=(mAG2_+tresh2)*0.5*(NN_+MM_) +...
	    (tresh2+mIG2_)*0.5*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/(1-lambda)<mIG2_)*RR_+...
        triangle(tresh2_,tresh2,tresh1_,tresh1,0,0)*PP_; %todo - check triangle mAL_kes sense

    vbar=Pww_*(vww) + Pwn_*vwn + Pnw_*vnw + lambda*conexp1_+(1-lambda)*conexp2_;
 
    
    V=vbar;
    conexp1=conexp1_;
    conexp2=conexp2_;

    Pnw=Pnw_;
    Pwn=Pwn_;
    Pww=Pww_;
    
    Pw0w0=0;
    Pw0n=0;
    Pw0w=0;

    Pnw0=0;
    Pww0=0;
        
    Pnn=0;
    workh=Pww+Pwn+Pw0w0+Pw0n+Pw0w+Pww0;
    workw=Pww+Pnw+Pw0w0+Pnw0+Pw0w+Pww0;

    check=round(Pww+Pwn+Pw0w0+Pw0n+Pw0w+Pww0+Pnw+Pnw0,2);
    if (round(Pww_ + Pwn_+ Pnw_,2)==1)...
    *(check==1)*(Pww<=1)...
    ~=1
        (Pww_ + Pwn_+ Pnw_)
        check
        Pww<=1
        %(conexp1<max(mA1_,mAG1_))
        %(conexp1>=min(mI1_,mIG1_)*workh)
        %(conexp2<max(mA2_,mAG2_))
        %(conexp2>=min(mI2_,mIG2_)*workw)
        %(conexp1<max(mA1_,mAG1_))*(conexp1>=min(mI1_,mIG1_)*workh)*(conexp2<max(mA2_,mAG2_))*(conexp2>=min(mI2_,mIG2_)*workw)
        return
    end
    %workh
    %workw

else
    Pwn=0;
    Pww=0;
    Pnw=0;
    Pww0=0;
    Pw0w=0;
    if vww>vwn && vww>vnw
        Pww=1;
        V=vww;
        workh=1;
        workw=1;
    elseif vwn>vww && vwn>vnw
        Pwn=1;
        V=vwn;
        workh=1;
        workw=0;
    elseif vnw>vww && vnw>vwn
        Pnw=1;
        V=vnw;
        workh=0;
        workw=1;
    elseif vnw==vww && vnw>vwn
        Pnw=0.5;
        Pww=0.5;
        V=vnw*Pnw+vww*Pww; 
        workh=0.5;
        workw=1;
    elseif vwn==vww && vwn>vnw
        Pwn=0.5;
        Pww=0.5;
        V=vwn*Pwn+vww*Pww; 
        workh=1;
        workw=0.5;
    elseif vnw==vwn && vnw>vww
        Pnw=0.5;
        Pwn=0.5;
        V=vnw*Pnw+vwn*Pwn; 
        workh=0.5;
        workw=0.5;
    end
   conexp1=0;
   conexp2=0;
   Pnw0=0;
   Pw0n=0;
   Pw0w0=0;
   Pnn=0;   

end

end
