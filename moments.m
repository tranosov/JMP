
function [moments_,other_,outcomes_]...
    =moments(p,EQS,PARREST, DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC, DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS,fast)

% VC NOT USED

global sizehs

%todo: only compute the necessary moments!

D=PARREST.('D');
%sigmal=PARREST.('sigmal');
sigmam=PARREST.('sigmam');
JLs=PARREST.('JLs');
mm=PARREST.('mm');
Jw=PARREST.('Jw');
%NC=PARREST.('NC');
XI=PARREST.('XI');
mA=PARREST.('mA');
mI=PARREST.('mI');
mAL=PARREST.('mAL');
mIL=PARREST.('mIL');
plocal=PARREST.('plocal');
typeic=PARREST.('typeic');
betah=PARREST.('betah');
AC=PARREST.('AC');
AS=PARREST.('AS');
HS=PARREST.('HS');
Jc=PARREST.('Jc');

params=PARREST.('params');
NLY=params{'NLY',:};

N=size(D);
I=N(1);
N=size(Jw);
T=N(2);
W=3;
sizehs=size(HS);

w1=EQS.('w1') ;
w2=EQS.('w2') ;
lei=EQS.('leis');

Yss=OUTS.('Ys') ;
lss=OUTS.('lss') ;
ucls=OUTS.('ucls') ;
uls=OUTS.('uls') ;
ucs=OUTS.('ucs') ;
uhs=OUTS.('uhs') ;
uxs=OUTS.('uxs') ;
xs=OUTS.('xs') ;

Yc=OUTC.('Yc') ;
xh=OUTC.('xh') ;
xw=OUTC.('xw') ;
lh=OUTC.('lh') ;
lw=OUTC.('lw') ;
ch=OUTC.('ch') ;
cw=OUTC.('cw') ;
uh=OUTC.('uh') ;
uw=OUTC.('uw') ;
uch=OUTC.('uch') ;
ucw=OUTC.('ucw') ;
ulh=OUTC.('ulh') ;
ulw=OUTC.('ulw') ;
ux=OUTC.('ux') ;
uho=OUTC.('uho') ;
ucl=OUTC.('ucl') ; % probably not nec here?

%vc=OUTC.('vc') ;- not used!

distS=zeros(I,T,I,W); % offered distance!
commuteS=zeros(I,T,I,W); %working
lsS=zeros(I,T,I,W);
worksS=zeros(I,T,I,W);
closeS=zeros(I,T,I,W);
YS=zeros(I,T,I,W);
L_low=zeros(I,T,I,W); % is the local location NOT an industry center?
matchS=zeros(I,T,I,W); 
uconsS=zeros(I,T,I,W);
uleisS=zeros(I,T,I,W);
uhousS=zeros(I,T,I,W);
uxS=zeros(I,T,I,W);
xS=zeros(I,T,I,W);
aS=zeros(I,T,I,W);
for j=1:I
    for t=1:T
        for i=1:I
            for w=2:W
                j_=j*(w==3) + i*(w==2);
                worksS(j,t,i,w)=(w>1);
                distS(j,t,i,w)=D(i,j); % offer
                L_low(j,t,i,w)=out5(@()matchdist(i,j,t,mA,mI,mAL,mIL,typeic,D,mm,JLs,betah));
                commuteS(j,t,i,w)=D(i,j_);
                lsS(j,t,i,w)=lss(t,j_,i);
                YS(j,t,i,w)=Yss(t,j_,i);
                closeS(j,t,i,w)=ucls(t,j_,i);
                %PHID*(FdistS(commuteS(j,t,i,w))-1);

          
                [~,~,~,~,~,low]=matchdist(i,j_,t,mA,mI,mAL,mIL,typeic,D,mm,JLs,betah); % of what you really took
                ic=1-low; % is j low (not an industry center)?
                matchS(j,t,i,w)=cexps(j,t,i)+ XI* lsS(j,t,i,w)*ic; %first part: is not really correct for w=2 vs w=3, but should work on the aggregate
                % i.e. the first part you can take as an expectation before
                % you see the match shocks!
                
                uconsS(j,t,i,w)=ucs(t,j_,i);
                %ceh*uu(cs(YS(j,t,i,w),p(i)));
                uleisS(j,t,i,w)=uls(t,j_,i); %uu(-alphah* lsS(j,t,i,w) - betah* commuteS(j,t,i,w)+1+tih);
                uhousS(j,t,i,w)=uhs(t,j_,i); %eta*uu(hdS(YS(j,t,i,w),p(i)));
                uxS(j,t,i,w)=uxs(t,j_,i);
                xS(j,t,i,w)=xs(t,j_,i);
                aS(j,t,i,w)=AS(i);
                wagehS(j,t,i,w)=w1(lsS(j,t,i,w),ic); 
            end
        end
    end
end


conshC=zeros(I,I,T,T,I,W,W); 
conswC=zeros(I,I,T,T,I,W,W);
uconshC=zeros(I,I,T,T,I,W,W); 
uconswC=zeros(I,I,T,T,I,W,W);
uhousC=zeros(I,I,T,T,I,W,W);
aC=zeros(I,I,T,T,I,W,W);
uhC_det=zeros(I,I,T,T,I,W,W); 
uwC_det=zeros(I,I,T,T,I,W,W);
uhC=zeros(I,I,T,T,I); 
uwC=zeros(I,I,T,T,I);
matchhC=zeros(I,I,T,T,I,W,W); 
matchwC=zeros(I,I,T,T,I,W,W);
lshC=zeros(I,I,T,T,I,W,W); 
lswC=zeros(I,I,T,T,I,W,W);
xhC=zeros(I,I,T,T,I,W,W); 
xwC=zeros(I,I,T,T,I,W,W);
commutehC=zeros(I,I,T,T,I,W,W); % misc
commutewC=zeros(I,I,T,T,I,W,W);
disthC=zeros(I,I,T,T,I,W,W); %offered
distwC=zeros(I,I,T,T,I,W,W);
leisurehC=zeros(I,I,T,T,I,W,W); 
leisurewC=zeros(I,I,T,T,I,W,W); 
leishC=zeros(I,I,T,T,I,W,W); 
leiswC=zeros(I,I,T,T,I,W,W);
closeC=zeros(I,I,T,T,I,W,W); 
publicC=zeros(I,I,T,T,I,W,W);
workshC=zeros(I,I,T,T,I,W,W); 
workswC=zeros(I,I,T,T,I,W,W); 
PwwC=zeros(I,I,T,T,I,W,W);
PwwC_1accept=zeros(I,I,T,T,I,W,W);
PwwC_0accept=zeros(I,I,T,T,I,W,W);
PwwC_2accept=zeros(I,I,T,T,I,W,W);
YC=zeros(I,I,T,T,I,W,W);
mind=zeros(I,I,T,T,I,W,W);
wagehC=zeros(I,I,T,T,I,W,W);
wagewC=zeros(I,I,T,T,I,W,W);
L_lowhC=zeros(I,I,T,T,I,W,W);
L_lowwC=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:3 
                       for ww=1:3
hi=wh>1;
wi=ww>1;
wwi=hi*(1-wi)+2*wi*(1-hi) + 3*wi*hi;
if wwi>0
    dh=D(i,jh)*(wh==3) + D(i,i)*(wh==2);
    dw=D(i,jw)*(ww==3) + D(i,i)*(ww==2);
    jh_=jh*(wh==3) + i*(wh==2|wh==1);
    jw_=jw*(ww==3) + i*(ww==2|ww==1);
    [~,~,~,~,~,low]=matchdist(i,jh_,th,mA,mI,mAL,mIL,typeic,D,mm,JLs,betah);
    ich=1-low;
    [~,~,~,~,~,low]=matchdist(i,jw_,tw,mA,mI,mAL,mIL,typeic,D,mm,JLs,betah);
    icw=1-low;
    lsh_=lh(th,tw,jh_,jw_,i,wwi);  
    lsw_=lw(th,tw,jh_,jw_,i,wwi);
   
    conshC(jh,jw,th,tw,i,wh,ww)=ch(th,tw,jh_,jw_,i,wwi);
    conswC(jh,jw,th,tw,i,wh,ww)=cw(th,tw,jh_,jw_,i,wwi);
    uconshC(jh,jw,th,tw,i,wh,ww)=uch(th,tw,jh_,jw_,i,wwi);
    uconswC(jh,jw,th,tw,i,wh,ww)=ucw(th,tw,jh_,jw_,i,wwi);
    lshC(jh,jw,th,tw,i,wh,ww)=lsh_; 
    lswC(jh,jw,th,tw,i,wh,ww)=lsw_;
    xhC(jh,jw,th,tw,i,wh,ww)=xh(th,tw,jh_,jw_,i,wwi);
    xwC(jh,jw,th,tw,i,wh,ww)=xw(th,tw,jh_,jw_,i,wwi);
    
    publicC(jh,jw,th,tw,i,wh,ww)=ux(th,tw,jh_,jw_,i,wwi);
    closeC(jh,jw,th,tw,i,wh,ww)=ucl(th,tw,jh_,jw_,i,wwi);
    uhC_det(jh,jw,th,tw,i,wh,ww)=uh(th,tw,jh_,jw_,i,wwi) ; % missing conexp
    uwC_det(jh,jw,th,tw,i,wh,ww)=uw(th,tw,jh_,jw_,i,wwi) ;
    disthC(jh,jw,th,tw,i,wh,ww)=D(i,jh); %to an offered job
    distwC(jh,jw,th,tw,i,wh,ww)=D(i,jw);
    commutehC(jh,jw,th,tw,i,wh,ww)=dh; %IF WORKING
    commutewC(jh,jw,th,tw,i,wh,ww)=dw;
    leisurehC(jh,jw,th,tw,i,wh,ww)=lei(lsh_,xhC(jh,jw,th,tw,i,wh,ww),dh);
    leisurewC(jh,jw,th,tw,i,wh,ww)=lei(lsw_,xwC(jh,jw,th,tw,i,wh,ww),dw);
    leishC(jh,jw,th,tw,i,wh,ww)=ulh(th,tw,jh_,jw_,i,wwi); %utility!
    leiswC(jh,jw,th,tw,i,wh,ww)=ulw(th,tw,jh_,jw_,i,wwi);
    workshC(jh,jw,th,tw,i,wh,ww)=hi;
    workswC(jh,jw,th,tw,i,wh,ww)=wi;
    PwwC(jh,jw,th,tw,i,wh,ww)=hi*wi;
    PwwC_1accept(jh,jw,th,tw,i,wh,ww)=(wh==3 & ww==2) | (wh==2 & ww==3);
    PwwC_2accept(jh,jw,th,tw,i,wh,ww)=(wh==3 & ww==3);
    PwwC_0accept(jh,jw,th,tw,i,wh,ww)=(wh==2 & ww==2);
    YC(jh,jw,th,tw,i,wh,ww)=Yc(th,tw,jh_,jw_,i,wwi);
    mind(jh,jw,th,tw,i,wh,ww)=min(dh*hi,dw*wi); %0 if not working
    wagehC(jh,jw,th,tw,i,wh,ww)=w1(lsh_,ich)*hi; %+99999*(1-hi); % 0 if not working
    wagewC(jh,jw,th,tw,i,wh,ww)=w2(lsw_,icw)*wi; %+99999*(1-wi); 
    L_lowhC(jh,jw,th,tw,i,wh,ww)=out5(@()matchdist(i,jh,th,mA,mI,mAL,mIL,typeic,D,mm,JLs,betah)); % is the local labor market 'low' for me?
    L_lowwC(jh,jw,th,tw,i,wh,ww)=out5(@()matchdist(i,jw,tw,mA,mI,mAL,mIL,typeic,D,mm,JLs,betah));

    uhousC(jh,jw,th,tw,i,wh,ww)=uho(th,tw,jh_,jw_,i,wwi);
    aC(jh,jw,th,tw,i,wh,ww)=AC(i);

    matchhC(jh,jw,th,tw,i,wh,ww)=  cexp1(jh,jw,th,tw,i) + XI*lsh_*ich ;
    matchwC(jh,jw,th,tw,i,wh,ww)=  cexp2(jh,jw,th,tw,i) + XI*lsw_*icw ;
end
                       end
                    end
        
        
		uhC(jh,jw,th,tw,i)= sum(sum(Pw(jh,jw,th,tw,i,:,:).*uhC_det(jh,jw,th,tw,i,:,:))) + cexp1(jh,jw,th,tw,i);
		uwC(jh,jw,th,tw,i)= sum(sum(Pw(jh,jw,th,tw,i,:,:).*uwC_det(jh,jw,th,tw,i,:,:))) + cexp2(jh,jw,th,tw,i);
                
                end
            end
        end
    end
end

%%

DSi=sumS(DS);
DCi=sumC(DC);

sDSi=DSi/sum(DSi);
sDCi=DCi/sum(DCi);
scity_=(DSi+DCi)/(sum(DSi)+sum(DCi));
scity=scity_(1);

%% rescale  so that distribution is what I want it to be:

%%
workhC_raw=sum(sumC(DC.*workshC))/sum(DCi);
workwC_raw=sum(sumC(DC.*workswC))/sum(DCi);

workS_raw=sum(sumS(DS.*worksS))/sum(DSi);
% hours - includes 0
hours0hC_location=sumC(DC.*lshC)./DCi;
hours0wC_location=sumC(DC.*lswC)./DCi;
hours0hC_suburb=(hours0hC_location(2)*DCi(2)+hours0hC_location(3)*DCi(3))/(DCi(2)+DCi(3));
hours0wC_suburb=(hours0wC_location(2)*DCi(2)+hours0wC_location(3)*DCi(3))/(DCi(2)+DCi(3));

%%
hours0hC_raw=sum(sumC(DC.*lshC))/sum(DCi); %0.25
hours0wC_raw=sum(sumC(DC.*lswC))/sum(DCi) ;%0.14


%%
comtimehC_raw=sum(sumC(DC.*commutehC.*workshC))./sum(sumC(DC));
comtimewC_raw=sum(sumC(DC.*commutewC.*workswC))./sum(sumC(DC));


commuteS_raw=sum(sumS(commuteS.*DS.*worksS))./sum(sumS(DS.*worksS));
commutehS_raw=sum(sumS(commuteS.*DSh.*worksS))./sum(sumS(DSh.*worksS));
commutewS_raw=sum(sumS(commuteS.*DSw.*worksS))./sum(sumS(DSw.*worksS));
commutehC_raw=sum(sumC(DC.*commutehC.*workshC))./sum(sumC(DC.*workshC));
commutewC_raw=sum(sumC(DC.*commutewC.*workswC))./sum(sumC(DC.*workswC));


% as usual - I need more commute differences and in general
%% utilities - frist from match
DC_agr=sum(sum(DC,7),6);
DC1_agr=sum(sum(DC1,7),6);
DC2_agr=sum(sum(DC2,7),6);
DS_agr=sum(DS,4);
DSh_agr=sum(DSh,4);
DSw_agr=sum(DSw,4);
matchhS_raw=sum(sumS(DSh.*matchS))./sum(sum(sum(DSh_agr)));
matchwS_raw=sum(sumS(DSw.*matchS))./sum(sum(sum(DSw_agr)));

matchhC_raw=sum(sumC(matchhC.*DC))./sum(DCi);
matchwC_raw=sum(sumC(matchwC.*DC))./sum(DCi);

%%
% hours - positive

hourshC_raw=sum(sumC(DC.*lshC.*workshC))./sum(sumC(DC.*workshC)); %0.26
hourswC_raw=sum(sumC(DC.*lswC.*workswC))./sum(sumC(DC.*workswC)) ;%0.18

%% if both work
hourshC_pww_raw=sum(sumC(DC.*lshC.*PwwC))./sum(sumC(DC.*PwwC));
hourswC_pww_raw=sum(sumC(DC.*lswC.*PwwC))./sum(sumC(DC.*PwwC));



%% if just husband works/ wife works
hourshC_h0_raw=sum(sumC(DC.*lshC.*workshC.*(1-workswC)))./sum(sumC(DC.*workshC.*(1-workswC)));
hourswC_0w_raw=sum(sumC(DC.*lswC.*workswC.*(1-workshC)))./sum(sumC(DC.*workswC.*(1-workshC)));
hourshC_h0_location=sumC(DC.*lshC.*workshC.*(1-workswC))./sumC(DC.*workshC.*(1-workswC));
hourswC_0w_location=sumC(DC.*lswC.*workswC.*(1-workshC))./sumC(DC.*workswC.*(1-workshC));

d_=sumC(DC.*workshC.*(1-workswC));
hourshC_h0_suburb=(hourshC_h0_location(2)*d_(2)+hourshC_h0_location(3)*d_(3))/(d_(2)+d_(3));
d_=sumC(DC.*workswC.*(1-workshC));
hourswC_0w_suburb=(hourswC_0w_location(2)*d_(2)+hourswC_0w_location(3)*d_(3))/(d_(2)+d_(3));

%% wages
lwagehC=log(wagehC);
lwagehC(isinf(lwagehC))=0;
lwagewC=log(wagewC);
lwagewC(isinf(lwagewC))=0;
lwagehC_raw=sum(sumC(DC.*lwagehC.*workshC))./sum(sumC(DC.*workshC));
lwagewC_raw=sum(sumC(DC.*lwagewC.*workswC))./sum(sumC(DC.*workswC)) ;
wagegap_hw=lwagehC_raw-lwagewC_raw;
wagegap_hw_within=sum(sumC(DC.*(lwagehC-lwagewC).*PwwC))./sum(sumC(DC.*PwwC));

lwagehS=log(wagehS);
lwagehS(isinf(lwagehS))=0;
lwagehS_raw=sum(sumS(DSh.*lwagehS))./sum(sumS(DSh)); 
lwagewS_raw=sum(sumS(DSw.*lwagehS))./sum(sumS(DSw)) ;

wagegap_actual=((lwagehS_raw*sum(sumS(DSh)) + lwagehC_raw*sum(sumC(DC.*workshC))))/(sum(sumS(DSh)) + sum(sumC(DC.*workshC))) ...
    - (lwagewS_raw*sum(sumS(DSw)) + lwagewC_raw*sum(sumC(DC.*workswC)))/(sum(sumS(DSw)) + sum(sumC(DC.*workswC)));

%%
% hours - singles

hours0S_raw=sum(sumS(lsS.*DS))./sum(sumS(DS));

%%


% consumption

conshC_raw=sum(sumC(DC.*conshC))/sum(DCi); 
conswC_raw=sum(sumC(DC.*conswC))/sum(DCi); %
 
uconshC_raw=sum(sumC(DC.*uconshC))/sum(DCi); 
uconswC_raw=sum(sumC(DC.*uconswC))/sum(DCi); 
%%
%conshS_location=sumS(DS.*conshS)./DSi;
%conswS_location=sumS(DS.*conswS)./DSi;
uconshS_raw=sum(sumS(DSh.*uconsS))/sum(sumS(DSh));
uconswS_raw=sum(sumS(DSw.*uconsS))/sum(sumS(DSh)); 
uleishS_raw=sum(sumS(DSh.*uleisS))/sum(sumS(DSh));
uleiswS_raw=sum(sumS(DSw.*uleisS))/sum(sumS(DSw)); 
%conshS_raw=sum(sumS(DSh.*conshS))/sum(DSi);
%conswS_raw=sum(sumS(DSw.*conswS))/sum(DSi); %


leisurehC_raw=sum(sumC(DC.*leisurehC))/sum(DCi); 
leisurewC_raw=sum(sumC(DC.*leisurewC))/sum(DCi); %


leishC_raw=sum(sumC(DC.*leishC))/sum(DCi); 
leiswC_raw=sum(sumC(DC.*leiswC))/sum(DCi); %

%%

xhC_raw=sum(sumC(DC.*xhC))/sum(DCi) ;
xwC_raw=sum(sumC(DC.*xwC))/sum(DCi) ;% a little low?

xhC_pww_raw=sum(sumC(DC.*xhC.*PwwC))/sum(sumC(DC.*PwwC)) ;
xwC_pww_raw=sum(sumC(DC.*xwC.*PwwC))/sum(sumC(DC.*PwwC)) ;
xhC_h0_raw=sum(sumC(DC.*xhC.*(1-workswC).*workshC))/sum(sumC(DC.*(1-workswC).*workshC)) ;
xwC_h0_raw=sum(sumC(DC.*xwC.*(1-workswC).*workshC))/sum(sumC(DC.*(1-workswC).*workshC)) ;
xhC_0w_raw=sum(sumC(DC.*xhC.*(1-workshC).*workswC))/sum(sumC(DC.*(1-workshC).*workswC)) ;
xwC_0w_raw=sum(sumC(DC.*xwC.*(1-workshC).*workswC))/sum(sumC(DC.*(1-workshC).*workswC)) ;


xwS_raw=sum(sumS(DSh.*xS))/sum(sumS(DSh));
xhS_raw=sum(sumS(DSw.*xS))/sum(sumS(DSw));
xS_raw=sum(sumS(DS.*xS))/sum(sumS(DS));
%%
publicwS_raw=sum(sumS(DSh.*uxS))/sum(sumS(DSh));
publichS_raw=sum(sumS(DSw.*uxS))/sum(sumS(DSw));

closeC_raw=sum(sumC(DC.*closeC))/sum(DCi);
publicC_raw=sum(sumC(DC.*publicC))/sum(DCi);


closehS_raw=sum(sumS(DSh.*closeS))/sum(sumS(DSh));


closewS_raw=sum(sumS(DSw.*closeS))/sum(sumS(DSw));



%%


uhS_life=2*sum(sum(sum(DSh_agr.*VS)))./sum(sum(sum(DSh_agr)));
uwS_life=2*sum(sum(sum(DSw_agr.*VS)))./sum(sum(sum(DSw_agr)));

uhS_raw=uhS_life/2;
uwS_raw=uwS_life/2;

VS_raw=sum(sum(sum(DS_agr.*VS)))./sum(sum(sum(DS_agr)));

uhC_raw=sum(sum(sum(sum(sum(DC_agr.*uhC)))))./sum(sum(sum(sum(sum(DC_agr)))));
uwC_raw=sum(sum(sum(sum(sum(DC_agr.*uwC)))))./sum(sum(sum(sum(sum(DC_agr)))));


uhC_life=sum(sum(sum(sum(sum(DC1_agr.*uhC)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uhC)))))./sum(sum(sum(sum(sum(DC2_agr)))));
uwC_life=sum(sum(sum(sum(sum(DC1_agr.*uwC)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uwC)))))./sum(sum(sum(sum(sum(DC2_agr)))));

wVMAR=uwC_life - uwS_life; 
hVMAR=uhC_life - uhS_life;


%%
 a=(DS.*worksS);
JSingle=sum(sum(a(:,:,:,3),4),3)+transpose(reshape(sum(sum(a(:,:,:,2),4),1),[T,I]));

a=DC.*workshC;
Jcoupleh=reshape(sum(sum(sum(sum(sum(a(:,:,:,:,:,3,:),7),6),5),2),4),[I,T])+...
    transpose(reshape(sum(sum(sum(sum(sum(a(:,:,:,:,:,2,:),7),6),1),2),4),[T,I]));

a=DC.*workswC; %TR: corrected recently!
Jcouplew=reshape(sum(sum(sum(sum(sum(a(:,:,:,:,:,:,3),7),6),5),1),3),[I,T])+...
    transpose(reshape(sum(sum(sum(sum(sum(a(:,:,:,:,:,:,2),7),6),2),1),3),[T,I]));

Jobs=JSingle+Jcoupleh+Jcouplew;

JLs_m=Jobs./repmat(sum(Jobs),I,1);
JLs_all_m=sum(Jobs,2)/sum(sum(Jobs));

%%
distoS_m=zeros(I,T,I,W); % singles in categories
distalljS_m=zeros(I,T,I,W); 
for j=1:I
    for t=1:T
        for i=1:I
            for w=1:W
                distoS_m(j,t,i,w)=Do(i,JLs_m(:,t),D); %indep of jh,jw
                distalljS_m(j,t,i,w)=Do(i,JLs_all_m,D);
            end
        end
    end
end


distoS_m_raw=sum(sumS(distoS_m.*DS))/sum(sumS(DS));
distohS_m_raw=sum(sumS(distoS_m.*DSh))/sum(sumS(DSh));
distalljS_m_raw=sum(sumS(distalljS_m.*DS))/sum(sumS(DS));

distohC_m=zeros(I,I,T,T,I,W,W);
distowC_m=zeros(I,I,T,T,I,W,W);
distalljC_m=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:W
                        for ww=1:W
distohC_m(jh,jw,th,tw,i,wh,ww)=  Do(i,JLs_m(:,th),D);
distowC_m(jh,jw,th,tw,i,wh,ww)=  Do(i,JLs_m(:,tw),D);
distalljC_m(jh,jw,th,tw,i,wh,ww)=  Do(i,JLs_all_m,D);
                        end
                    end
                end
            end
        end
    end
end   


distohC_m_raw=sum(sumC(distohC_m.*DC))./sum(DCi);
distowC_m_raw=sum(sumC(distowC_m.*DC))./sum(DCi);
distalljC_m_raw=sum(sumC(distalljC_m.*DC))./sum(DCi);

disthC_raw=sum(sumC(DC.*disthC))/sum(DCi); % OFFERS
distwC_raw=sum(sumC(DC.*distwC))/sum(DCi);

%% standard deviation of do and dj

do_mean=(sum(sumC(DC.*(distohC_m))) + sum(sumC(DC.*(distowC_m))) + sum(sumS(DSh.*(distoS_m))) + sum(sumS(DSw.*(distoS_m))))/(sum(sumC(DC))*2 + sum(sumS(DS)));
do_sd=(sum(sumC(DC.*(distohC_m -do_mean).^2)) + sum(sumC(DC.*(distowC_m -do_mean).^2)) + sum(sumS(DSh.*(distoS_m -do_mean).^2)) + sum(sumS(DSw.*(distoS_m -do_mean).^2)))/(sum(sumC(DC))*2 + sum(sumS(DS)));

dj_mean=(sum(sumC(DC.*(distalljC_m))) + sum(sumC(DC.*(distalljC_m))) + sum(sumS(DSh.*(distalljS_m))) + sum(sumS(DSw.*(distalljS_m))))/(sum(sumC(DC))*2 + sum(sumS(DS)));
dj_sd=(sum(sumC(DC.*(distalljC_m -dj_mean).^2)) + sum(sumC(DC.*(distalljC_m -dj_mean).^2)) + sum(sumS(DSh.*(distalljS_m -dj_mean).^2)) + sum(sumS(DSw.*(distalljS_m -dj_mean).^2)))/(sum(sumC(DC))*2 + sum(sumS(DS)));


%% and of d itself!

d_mean=(sum(sumC(DC.*(commutehC.*workshC))) + sum(sumC(DC.*(commutewC.*workswC))) + sum(sumS(DSh.*(commuteS))) + sum(sumS(DSw.*(commuteS))))/(sum(sumC(DC))*2 + sum(sumS(DS)));
d_sd=(sum(sumC(DC.*(commutehC-d_mean).*workshC.^2)) + sum(sumC(DC.*(commutewC-d_mean).*workswC.^2)) + sum(sumS(DSh.*(commuteS-d_mean).^2)) + sum(sumS(DSw.*(commuteS-d_mean).^2)))/(sum(sumC(DC))*2 + sum(sumS(DS)));



%% abs hdo-wdo
abs_hdo_wdo=sum(sumC(abs(distowC_m-distohC_m).*DC))./sum(DCi);
mean_hdo_wdo=sum(sumC((distowC_m-distohC_m).*DC))./sum(DCi);
sd_hdo_wdo=(sum(sumC(DC.*((distowC_m-distohC_m-mean_hdo_wdo).^2))))/(sum(sumC(DC)));


%% jobjob distances
jobjob_prep_m=zeros(I,T+1);
for i=1:I
    for t=1:T
        jobjob_prep_m(i,t)=Do(i,JLs_m(:,t),D); %indep of jh,jw, based on all jobs, just by type
    end
    jobjob_prep_m(i,T+1)=Do(i,JLs_all_m,D);
end

jobjob_all_m=sum(jobjob_prep_m(:,T+1).*JLs_all_m)/sum(JLs_all_m) ; % data 10.7
jobjob_within_m=sum(sum(jobjob_prep_m(:,1:T).*JLs_m))/sum(sum(JLs_m)); % data 10.2


JLsh_m=Jcoupleh./repmat(sum(Jcoupleh),I,1);
jobjob_hwithin_m=sum(sum(jobjob_prep_m(:,1:T).*JLsh_m))/sum(sum(JLsh_m)); % need to replace JLs_m with dist of husband jobs?

% hw potential job distances - since matching random
jobjob_inds_m=zeros(T,T);
for t1=1:T
    for t2=1:T
        jobjob_inds_m(t1,t2)=sum(jobjob_prep_m(:,t1).*JLs_m(:,t2))/sum(JLs_m(:,t2));
    end
end
jobjob_hw_m=sum(sum(jobjob_inds_m.*Jc))/sum(sum(Jc));

% hw job distance? (not offer, actual job!) 
jobjob_hw_toact=zeros(I,I,T,T,I,W,W);
jobjob_hw_tooff=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:W
                        for ww=1:W
jh_=jh*(wh==3) + i*(wh==2|wh==1);
jw_=jw*(ww==3) + i*(ww==2|ww==1);

jobjob_hw_toact(jh,jw,th,tw,i,wh,ww)=  D(jh_,jw_);                            
jobjob_hw_tooff(jh,jw,th,tw,i,wh,ww)=  D(jh,jw);
                        end
                    end
                end
            end
        end
    end
end   
jobjob_hwact=sum(sumC(jobjob_hw_toact.*DC.*PwwC))./sum(sumC(DC.*PwwC)) ;
% Housing demand
%%

% share of income spent on housing
pHSi=sumS(HSingle).*p;
pHCi=sumC(HCouple).*p;
sHS=sum(pHSi)/sum(sumS(YS.*DS));
sHC=sum(pHCi)/sum(sumC(YC.*DC));
sH=(sum(pHSi)+sum(pHCi))/(sum(sumS(YS.*DS))+sum(sumC(YC.*DC)));

YncomeC=sum(sumC(YC.*DC))./sum(DCi);
YncomeS=sum(sumS(YS.*DS))./sum(DSi);
hC=sum(sumC(HCouple))/sum(DCi);


NLY_=(NLY*0.5*sum(sumS(DS))+ NLY*sum(sumC(DC)))/(sum(sumS(YS.*DS))+sum(sumC(YC.*DC)));
%%
% distance to opportunities vs working vs commuting
DC_long=reshape(DC.*(DC>=0),I*I*T*T*I*W*W,1);
DC_worksboth_long=reshape(DC.*PwwC.*(DC>=0),I*I*T*T*I*W*W,1);
worksh_long=reshape(workshC,I*I*T*T*I*W*W,1);
worksw_long=reshape(workswC,I*I*T*T*I*W*W,1);
lsh_long=8760*reshape(lshC,I*I*T*T*I*W*W,1); % good for when both work as well I think
lsw_long=8760*reshape(lswC,I*I*T*T*I*W*W,1);
commsh_long=reshape(commutehC,I*I*T*T*I*W*W,1); % needs to be when both work!
commsw_long=reshape(commutewC,I*I*T*T*I*W*W,1);
%dosh_long=reshape(distohC,I*I*T*T*I*W*W,1);
%dosw_long=reshape(distowC,I*I*T*T*I*W*W,1);
dosh_m_long=reshape(distohC_m,I*I*T*T*I*W*W,1);
dosw_m_long=reshape(distowC_m,I*I*T*T*I*W*W,1);

disth_long=reshape(disthC,I*I*T*T*I*W*W,1);
distw_long=reshape(distwC,I*I*T*T*I*W*W,1);
housh_long=8760*reshape(xhC,I*I*T*T*I*W*W,1);
housw_long=8760*reshape(xwC,I*I*T*T*I*W*W,1);
wageh_long=reshape(wagehC,I*I*T*T*I*W*W,1);
wagew_long=reshape(wagewC,I*I*T*T*I*W*W,1);
wageh_long(isinf(wageh_long))=0;
wagew_long(isinf(wagew_long))=0;

% with couple fes? how? would it change?
% TODO: so far I think these regressions are nonsense. distoh is by
% definition CONSTANT within couple, because I have both kinds of jobs
% being the same
location=zeros(I,I,T,T,I,W,W);
typeh=zeros(I,I,T,T,I,W,W);
typew=zeros(I,I,T,T,I,W,W);
loh=zeros(I,I,T,T,I,W,W);
low=zeros(I,I,T,T,I,W,W);
accepth=zeros(I,I,T,T,I,W,W);
acceptw=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:W
                        for ww=1:W
                location(jh,jw,th,tw,i,wh,ww)=i;
                typeh(jh,jw,th,tw,i,wh,ww)=th;
                typew(jh,jw,th,tw,i,wh,ww)=tw;
                loh(jh,jw,th,tw,i,wh,ww)=jh;
                low(jh,jw,th,tw,i,wh,ww)=jw;
                accepth(jh,jw,th,tw,i,wh,ww)=wh;
                acceptw(jh,jw,th,tw,i,wh,ww)=ww;
                        end
                    end
                end
            end
        end
    end
end
location_long=reshape(location,I*I*T*T*I*W*W,1);
typeh_long=reshape(typeh,I*I*T*T*I*W*W,1);
typew_long=reshape(typew,I*I*T*T*I*W*W,1);
loh_long=reshape(loh,I*I*T*T*I*W*W,1);
low_long=reshape(low,I*I*T*T*I*W*W,1);
accepth_long=reshape(accepth,I*I*T*T*I*W*W,1);
acceptw_long=reshape(acceptw,I*I*T*T*I*W*W,1);
x=1;
groups=zeros(I*I*T*T*I*W*W,1); % does it matter if I have only so many couples each with a weight?
for j=1:I*I*T*T*I*W*W
    groups(j)=x;
    x=x+1;
end

groups=[groups;groups];
weights=[DC_long;DC_long];
weights_bothwork=[DC_worksboth_long;DC_worksboth_long];

tbl=mat2dataset([[dosh_m_long;dosw_m_long],groups,[lsh_long;lsw_long],[commsh_long;commsw_long], ...
    [zeros(I*I*T*T*I*W*W,1);ones(I*I*T*T*I*W*W,1)],[zeros(I*I*T*T*I*W*W,1);dosw_m_long],weights,weights_bothwork,...
    [disth_long;distw_long],[zeros(I*I*T*T*I*W*W,1);distw_long],[location_long;location_long],...
    [typeh_long;typew_long],[typeh_long;typeh_long],[typew_long;typew_long],[loh_long;loh_long],[low_long;low_long],...
    [accepth_long;accepth_long],[acceptw_long;acceptw_long],[worksh_long;worksw_long],...
    [housh_long;housw_long],[wageh_long; wagew_long]]);


tbl.Properties.VarNames = {'dos_m','couple','ls','comm','woman','womandos_m',...
    'weights','weights_bothwork','dj','womandj','location','type','typeh','typew','loh','low',...
    'accepth','acceptw','works','hous','wage'};
tbl.couple=categorical(tbl.couple);

typedif=[[typeh_long==1] - [typew_long==1]];
lsdif=[lsh_long- lsw_long];
womdif=[zeros(I*I*T*T*I*W*W,1) - ones(I*I*T*T*I*W*W,1)];
dodif=[dosh_m_long - dosw_m_long];
wdodif=[zeros(I*I*T*T*I*W*W,1) - dosw_m_long];
comdif=[commsh_long - commsw_long];
workdif=[worksh_long - worksw_long];
hwkdif=[housh_long - housw_long];
wagedif=[wageh_long - wagew_long];

X=[typedif,womdif,dodif,wdodif];
Xplus=X(DC_worksboth_long>0,:);
X=X(DC_long>0,:);


WW=repmat(DC_long(DC_long>0),1,size(X,1)).*eye(size(X,1));
WWplus=repmat(DC_worksboth_long(DC_worksboth_long>0),1,size(Xplus,1)).*eye(size(Xplus,1));

betas=zeros(6,2);


try chol((X'*WW*X));
    y=lsdif(DC_long>0,:);
    b=((X'*WW*X)\eye(size(X'*WW*X)))*(X'*WW*y);
    betas(1,:)=[b(end-1),b(end)];
    
    y=workdif(DC_long>0,:);
    b=((X'*WW*X)\eye(size(X'*WW*X)))*(X'*WW*y);
    betas(3,:)=[b(end-1),b(end)];

    y=hwkdif(DC_long>0,:);
    b=((X'*WW*X)\eye(size(X'*WW*X)))*(X'*WW*y);
    betas(5,:)=[b(end-1),b(end)];
catch ME;
    betas(1,:)=10^12;
    betas(3,:)=10^12;
    betas(5,:)=10^12;
    disp('Matrix (X*WW*X) is not symmetric positive definite')
end


try chol((Xplus'*WWplus*Xplus));
   y=lsdif(DC_worksboth_long>0,:);
    b=((Xplus'*WWplus*Xplus)\eye(size(Xplus'*WWplus*Xplus)))*(Xplus'*WWplus*y);
    betas(2,:)=[b(end-1),b(end)];

    y=comdif(DC_worksboth_long>0,:);
    b=((Xplus'*WWplus*Xplus)\eye(size(Xplus'*WWplus*Xplus)))*(Xplus'*WWplus*y);
    betas(4,:)=[b(end-1),b(end)];

    y=wagedif(DC_worksboth_long>0,:);
    b=((Xplus'*WWplus*Xplus)\eye(size(Xplus'*WWplus*Xplus)))*(Xplus'*WWplus*y);
    betas(6,:)=[b(end-1),b(end)];
    
catch ME;
    disp('Matrix (X*WW*X)_plus is not symmetric positive definite')
    betas(2,:)=10^12;
    betas(4,:)=10^12;
    betas(6,:)=10^12;
end






%{
groupdummies=[eye(I*I*T*T*I*W*W);eye(I*I*T*T*I*W*W)];
typedummies=[[typeh_long==1];[typew_long==1]] % should be just 1 dummy? they are colinear?

X=[groupdummies(:,2:end), typedummies,[zeros(I*I*T*T*I*W*W,1);ones(I*I*T*T*I*W*W,1)],[dosh_m_long;dosw_m_long], [zeros(I*I*T*T*I*W*W,1);dosw_m_long]];
X=X(weights>0,:);
y=[lsh_long;lsw_long];
y=y(weights>0)
WW=repmat(weights(weights>0),1,size(X,1)).*eye(size(X,1));
XX=((X'*WW*X)\eye(size(X'*WW*X)))*(X'*WW*y)
%}
%% todo: rewrite for speed
%{
tic
betas_dos_m=zeros(6,4);
betas_womandos_m=zeros(6,4);
lm = fitlm(tbl,' ls ~ 1+ dos_m +woman+ womandos_m + couple+type','Weights',weights); % a bit slow
betas_womandos_m(1,:)=table2array(lm.Coefficients('womandos_m',:));
betas_dos_m(1,:)=table2array(lm.Coefficients('dos_m',:));
w = warning('query','last');
id = w.identifier;
warning('off',id)

% this one is 'too strong'
lm = fitlm(tbl,' ls ~ 1+ dos_m +woman+ womandos_m + couple+type','Weights',weights_bothwork);
betas_womandos_m(2,:)=table2array(lm.Coefficients('womandos_m',:));
betas_dos_m(2,:)=table2array(lm.Coefficients('dos_m',:));
lm = fitlm(tbl,' works ~ 1+ dos_m +woman+ womandos_m + couple+type','Weights',weights);
betas_womandos_m(3,:)=table2array(lm.Coefficients('womandos_m',:));
betas_dos_m(3,:)=table2array(lm.Coefficients('dos_m',:));
% kind of (a bit too big?).

lm = fitlm(tbl,' comm ~ 1+ dos_m +woman+ womandos_m + couple+type','Weights',weights_bothwork);
betas_womandos_m(4,:)=table2array(lm.Coefficients('womandos_m',:));
betas_dos_m(4,:)=table2array(lm.Coefficients('dos_m',:));

lm = fitlm(tbl,' hous ~ 1+ dos_m +woman+womandos_m+ couple+type','Weights',weights);
betas_womandos_m(5,:)=table2array(lm.Coefficients('womandos_m',:)); %hours (and it is causd by intensive margin) do not seem MORE sensitive for women?? CHECK
betas_dos_m(5,:)=table2array(lm.Coefficients('dos_m',:)); % only through extensive margin! with both working - commuting - less time goes the other way!


lm = fitlm(tbl,' wage ~ 1+ dos_m +woman+womandos_m+ couple+type','Weights',weights_bothwork);
betas_womandos_m(6,:)=table2array(lm.Coefficients('womandos_m',:))*100; %hours (and it is causd by intensive margin) do not seem MORE sensitive for women?? CHECK
betas_dos_m(6,:)=table2array(lm.Coefficients('dos_m',:));
toc
%}
%%
uhousC_raw=sum(sumC(uhousC.*DC))./sum(DCi);
uhoushS_raw=sum(sumS(DSh.*uhousS))./sum(sumS(DSh));
uhouswS_raw=sum(sumS(DSw.*uhousS))./sum(sumS(DSw));

%%
p_gradient_simple=log(p(1))-log(((p(2)+p(3))/2)); %city vs sub
dj=zeros(I,1);
for i=1:I
    dj(i)=Do(i,JLs_all_m,D);
end

employees=reshape(sum(sum(sum(sum(sum(sum(DC.*workshC)),4),5),6),7),[1,2])+...
    reshape(sum(sum(sum(sum(sum(sum(DC.*workswC)),3),5),6),7),[1,2])+...
    sum(sum(sum(DS.*worksS),3),4);

semployees=(reshape(sum(sum(sum(sum(sum(sum(DC.*workswC)),3),5),6),7),[1,2])+...
    sum(sum(sum(DSw.*worksS),3),4))./employees;

for i=1:I
    for t=1:T
        djj(i,t)=Do(i,JLs_m(:,t),D);
        prices(i,t)=log(p(i));
        sharew(i,t)=semployees(t);

    end
end

p_long=reshape(log(p),I,1);
dj_long=reshape(dj,I,1);

X=[ones(I,1),dj_long];
b=((X'*X)\eye(size(X'*X)))*(X'*p_long);
pg(1,:)=b(end);

%tbl=mat2dataset([p_long,dj_long]);
%tbl.Properties.VarNames = {'p','dj'};
%lm = fitlm(tbl,'  p ~ 1+ dj');
%pg(1,:)=table2array(lm.Coefficients('dj',:)); % log-point per mile

pp_long=reshape(prices,I*T,1);
djj_long=reshape(djj,I*T,1);
sharew_long=reshape(sharew,I*T,1);

X=[ones(I*T,1),djj_long,sharew_long,sharew_long.*djj_long];
b=((X'*X)\eye(size(X'*X)))*(X'*pp_long);
pg(2,:)=b(end);

%tbl=mat2dataset([pp_long,djj_long,sharew_long,sharew_long.*djj_long]);
%tbl.Properties.VarNames = {'p','dj','sw','sw_djj'};
%lm = fitlm(tbl,'  p ~ 1+ dj + sw +sw_djj');
%pg(2,:)=table2array(lm.Coefficients('sw_djj',:)); % log-point per mile

%%
distohC=zeros(I,I,T,T,I,W,W);
distowC=zeros(I,I,T,T,I,W,W);
distalljC=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:W
                        for ww=1:W
distohC(jh,jw,th,tw,i,wh,ww)=  Do(i,JLs(:,th),D);
distowC(jh,jw,th,tw,i,wh,ww)=  Do(i,JLs(:,tw),D);
                        end
                    end
                end
            end
        end
    end
end


distohC_raw=sum(sumC(distohC.*DC))./sum(DCi);
distowC_raw=sum(sumC(distowC.*DC))./sum(DCi);
%%

swnmarried=1 - (exp((wVMAR)/sigmam)/(1+exp((wVMAR)/sigmam)));
shnmarried=1 - (exp((hVMAR)/sigmam)/(1+exp((hVMAR)/sigmam)));
%%

Keys = {'scity_dif', 'wlfp_dif', 'hlfp','whours_pww_dif',...
    'scommiles','hdj','shcommiles_dif','swcommiles_difw',...
    'hrs0_dif_citysub',...
    'betacom' ,'betacom_w','betalfp' ,'betalfp_w','betahrs0','betahrs0_w','betahrs','betahrs_w','betahwk','betahwk_w',...
    'hdo','wdo_dif','sdo_dif','commiles_sd','do_sd','dj_sd','betalwg','betalwg_w',...
    'shours','hhours_pww','hhours_pwn_dif','whours_pnw_difwpww',...
    'whwk_pww','whwk_pwn_dif','whwk_pnw_dif','hhwk_pww_dif','hhwk_pwn_difhpww','hhwk_pnw_difhpww',...
    'd_jobjob_all_psid','d_jobjob_within_psid','d_jobjob_hwithin_psid','d_jobjob_hw','d_jobjob_hw_actual',...
    'abs_hdo_wdo','shouseexp','p_gradient','wagegap_hw_withn','shwk',...
    'snmarried','scity','sjobscity','sd_hdo_wdo','p_gradient_simple','NLY_'};
Keyso =  {'shnmarried','swnmarried','sdj_dif','sshouseexp','cshouseexp','wagegap_actual','wagegap_hw',...
    'p_gradient_swdjobs',...
    'hVC_hVS' , 'wVC_wVS','hL','wL','hl','wl','hx','wx','hcom0','wcom0','hcom','wcom','hcons','wcons',...
    'uconshC_raw','uconswC_raw','uconshS_raw','uconswS_raw',...
    'uleishC_raw','uleiswC_raw','uleishS_raw','uleiswS_raw',...
    'matchhC_raw','matchwC_raw','matchhS_raw','matchwS_raw',...
    'uhousC_raw','uhoushS_raw','uhouswS_raw',...
    'publicC_raw','closeC_raw','closehS_raw','closewS_raw','publichS_raw','publicwS_raw',...
    'hdo_of','wdo_of_dif','wdo_act_dif','Yc','HC','p1','p2','p3'};

model_m=[sDSi(1)-sDCi(1),workwC_raw-workhC_raw, workhC_raw, (hourswC_pww_raw-hourshC_pww_raw)*24*365,...
        commuteS_raw,distalljC_m_raw, commuteS_raw- commutehC_raw,commuteS_raw- commutewC_raw,...
        (hours0hC_suburb-hours0wC_suburb-(hours0hC_location(1)-hours0wC_location(1)))*24*365,...
        betas(4,1),betas(4,2), betas(3,1),betas(3,2),betas(1,1),betas(1,2),betas(2,1),betas(2,2),betas(5,1),betas(5,2),...
        distohC_m_raw,distowC_m_raw-distohC_m_raw,distohS_m_raw-distohC_m_raw, d_sd,do_sd, dj_sd, betas(6,1),betas(6,2),...
        hours0S_raw*24*365,hourshC_pww_raw*24*365,(hourshC_h0_raw-hourshC_pww_raw)*24*365,(hourswC_0w_raw-hourswC_pww_raw)*24*365,...
        xwC_pww_raw*24*365, ( xwC_h0_raw - xwC_pww_raw)*24*365, (xwC_0w_raw - xwC_pww_raw)*24*365, (xhC_pww_raw-xwC_pww_raw)*24*365,( xhC_h0_raw - xhC_pww_raw)*24*365, (xhC_0w_raw - xhC_pww_raw)*24*365,...
        jobjob_all_m, jobjob_within_m,jobjob_hwithin_m, jobjob_hw_m, jobjob_hwact,...
        abs_hdo_wdo,sH,  pg(1,1),-wagegap_hw_within,xS_raw*24*365,...
        shnmarried,scity,JLs_all_m(1,1),sd_hdo_wdo,p_gradient_simple,NLY_];
model_m_add=[shnmarried,swnmarried,distalljS_m_raw-distalljC_m_raw, sHS, sHC,wagegap_actual,wagegap_hw , pg(2,1) ];
other_=[hVMAR ,wVMAR, leisurehC_raw,leisurewC_raw,hourshC_raw,hourswC_raw,xhC_raw,xwC_raw,comtimehC_raw*betah,comtimewC_raw*betah,commutehC_raw*betah,commutewC_raw*betah,conshC_raw,conswC_raw,...
   uconshC_raw,uconswC_raw,uconshS_raw,uconswS_raw,...
    leishC_raw,leiswC_raw,uleishS_raw,uleiswS_raw,...
    matchhC_raw,matchwC_raw,matchhS_raw,matchwS_raw,...
    uhousC_raw,uhoushS_raw,uhouswS_raw,...
    publicC_raw,closeC_raw,closehS_raw,closewS_raw,publichS_raw,publicwS_raw,...
    distowC_raw,distowC_raw-distohC_raw,distwC_raw-disthC_raw,YncomeC, hC,p];

if I==4
    Keyso{end+1}='p4';
end
moments_= array2table([model_m]', 'RowNames',Keys);
other_= array2table([model_m_add,other_], 'VariableNames',Keyso);

aC_raw=sum(sumC(aC.*DC))./sum(DCi);
ahS_raw=sum(sumS(DSh.*aS))./sum(sumS(DSh));
awS_raw=sum(sumS(DSw.*aS))./sum(sumS(DSw));
outcomes_=[workhC_raw,workwC_raw,...
    sDCi(1),sDSi(1),...
    commutehC_raw,commutewC_raw,commuteS_raw,...
    distohC_m_raw,distowC_m_raw,distoS_m_raw,...
    uhC_raw,uwC_raw,uhS_raw,uwS_raw, VS_raw,...
    uconshC_raw,uconswC_raw,uconshS_raw,uconswS_raw,...
    leishC_raw,leiswC_raw,uleishS_raw,uleiswS_raw,...
    matchhC_raw,matchwC_raw,matchhS_raw,matchwS_raw,...
    uhousC_raw,uhoushS_raw,uhouswS_raw,...
    publicC_raw,closeC_raw,closehS_raw,closewS_raw,publichS_raw,publicwS_raw,...
    aC_raw,ahS_raw,awS_raw,(shnmarried+swnmarried)/2];
end