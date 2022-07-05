%{
Couples and singles location

Demand functions for couples.

Job offer in j X j X type t X t.
%}


function [DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC]=DCouple(p,alloutput,lambda,EQS,PARREST) 
%global D  AS  AC w1 w2 crra us ch cw cs hdS  Ys Yc lss_G lsah_G lsaw_G lsb_G alphaw betaw alphah betah tih tiw mu JLs Jw Jm HS Jc NC NS NSh NSw sigmaw sigmal mA mI kappa PHI mA mI mAL mIL
%global typeic OUTC EQS JLs D AC betah WARNINGS mm sigmal NC mu
global WARNINGS 

D=PARREST.('D');
sigmal=PARREST.('sigmal');
JLs=PARREST.('JLs');
mu=PARREST.('mu');
NC=PARREST.('NC');


I=size(D,1);
T=size(JLs,2);
W=3;

WARNINGS=0; % reset
OUTC=labors(p,alloutput,lambda,EQS,PARREST);
%OUTC_=OUTC;
%OUTC.vc

if WARNINGS>0
    DC=999*10^6;
    DC1=999;
    DC2=999;
    HCouple=999;
    Pw=999;
    hC=999;
    VC=999;
    cexp1=999;
    cexp2=999;
    OUTC=999;
    return
end
%{
vc=OUTC.('vc') ;
uh=OUTC.('uh');
uw=OUTC.('uw');
ulh=OUTC.('ulh');
ulw=OUTC.('ulw');
uch=OUTC.('uch');
ucw=OUTC.('ucw');
ich_=OUTC.('ich_');
icw_=OUTC.('icw_');
uXih=OUTC.('uXih');
uXiw=OUTC.('uXiw');
lh=OUTC.('lh');
lw=OUTC.('lw');
xh=OUTC.('xh');
xw=OUTC.('xw');
ux=OUTC.('ux');
ucl=OUTC.('ucl');
uho=OUTC.('uho');
% no by offers - as if all take offers
uh(1,2,2,3,2,3);
uh(1,2,2,3,3,3);
uh(1,2,2,3,2,3)-uh(1,2,2,3,3,3) % somehow - husbands utility is MORE negatively affected by couple settling in husband suburb. 
uw(1,2,2,3,2,3);
uw(1,2,2,3,3,3);
uw(1,2,2,3,2,3)-uw(1,2,2,3,3,3)

ulh(1,2,2,3,2,3)-ulh(1,2,2,3,3,3)
ulw(1,2,2,3,2,3)-ulw(1,2,2,3,3,3)
uXih(1,2,2,3,2,3)-uXih(1,2,2,3,3,3) 
uXiw(1,2,2,3,2,3)-uXiw(1,2,2,3,3,3) 
uch(1,2,2,3,2,3)-uch(1,2,2,3,3,3) 
ucw(1,2,2,3,2,3)-ucw(1,2,2,3,3,3) 
ux(1,2,2,3,2,3)-ux(1,2,2,3,3,3)
ucl(1,2,2,3,2,3)-ucl(1,2,2,3,3,3)
uho(1,2,2,3,2,3)-uho(1,2,2,3,3,3)
lh(1,2,2,3,2,3)-lh(1,2,2,3,3,3)
lw(1,2,2,3,2,3)-lw(1,2,2,3,3,3) 
%xh(1,2,2,3,2,3)-xh(1,2,2,3,3,3)
%xw(1,2,2,3,2,3)-xw(1,2,2,3,3,3) 


%ich_(1,2,2,3,2,3)-ich_(1,2,2,3,3,3) 
%icw_(1,2,2,3,2,3)-icw_(1,2,2,3,3,3)
lh(1,2,2,3,2,3)-lh(1,2,2,3,3,3)
lw(1,2,2,3,2,3)-lw(1,2,2,3,3,3) 

vc(1,2,2,3,2,3);
vc(1,2,2,3,3,3);
vc(1,2,2,3,2,3)-vc(1,2,2,3,3,3)

%vc(1,2,2,2,2,3);
%vc(1,2,3,3,3,3);
%vc(1,2,2,2,2,3)-vc(1,2,3,3,3,3) % this is driven by ich, byt it is way smaller, would never be enough.


EQS.('lambda')=0.51;
V=EQS.('V_lambda');
V=@(workh,workw,dh,dw,a,ch,cw,h,lsh,lsw,xh,xw,q,ich,icw) V(workh,workw,dh,dw,a,ch,cw,h,lsh,lsw,xh,xw,q,ich,icw,EQS.('lambda'));
hdC=EQS.('hdC');
lsh=EQS.('lsh') ;
lsw=EQS.('lsw');
ch_fun=EQS.('ch') ;
cw_fun=EQS.('cw');
xh_fun=EQS.('xh_fun');
xw_fun=EQS.('xw_fun');

vc=zeros(T,T,I,I,I,3); 
inputs=EQS.('inputs'); % last outputs
for i=1:I
    for jh=1:I
        for jw=1:I 
            for th=1:T
                for tw=1:T
                pp=p(i);
                dh=D(i,jh);
                dw=D(i,jw);
                [~,~,~,~,~,low]=matchdist(i,jh,th,0,0,0,0,typeic,D);
                ich=1-low;
                [~,~,~,~,~,low]=matchdist(i,jw,tw,0,0,0,0,typeic,D);
                icw=1-low;
                output1=inputs(th,tw,jh,jw,i,1,:);
                output2= inputs(th,tw,jh,jw,i,2,:);  
                output3=inputs(th,tw,jh,jw,i,3,:);
                vc(th,tw,jh,jw,i,1)=V(1,0,dh,0,AC(i),ch_fun(output1(1)),cw_fun(output1(1)),hdC(output1(1),pp),lsh(output1(1),output1(2),output1(3),dh,0),0,...
                   xh_fun(output1(1),output1(2),output1(3),dh,0),xw_fun(output1(1),output1(2),output1(3),dh,0),0,ich,0); 
                vc(th,tw,jh,jw,i,2)=V(0,1,0,dw,AC(i),ch_fun(output2(1)),cw_fun(output2(1)),hdC(output2(1),pp),0,lsw(output2(1),output2(2),output2(3),0,dw),...
                   xh_fun(output2(1),output2(2),output2(3),0,dw),xw_fun(output2(1),output2(2),output2(3),0,dw),0,0,icw); 
                vc(th,tw,jh,jw,i,3)=V(1,1,dh,dw,AC(i),ch_fun(output3(1)),cw_fun(output3(1)),hdC(output3(1),pp),lsh(output3(1),output3(2),output3(3),dh,dw),lsw(output3(1),output3(2),output3(3),dh,dw),...
                   xh_fun(output3(1),output3(2),output3(3),dh,dw),xw_fun(output3(1),output3(2),output3(3),dh,dw),0,ich,icw); 
               
                
                end
            end
        end
    end
end
OUTC.('vc') = vc;

vc(1,2,2,3,2,3);
vc(1,2,2,3,3,3);
vc(1,2,2,3,2,3)-vc(1,2,2,3,3,3) % bigger difference (i=3 BETTER more if husband valued more)


%vc(1,2,2,2,2,3);
%vc(1,2,3,3,3,3);
%vc(1,2,2,2,2,3)-vc(1,2,3,3,3,3); % this is driven by ich, byt it is way smaller, would never be enough - AND IT GETS SMALLER WITH LAMBDA higher???

%}
VC=zeros(I,I,T,T,I); %overall couple value
Pw=zeros(I,I,T,T,I,W,W);
VC1=zeros(I,I,T,T,I);
VC2=zeros(I,I,T,T,I);
cexp1=zeros(I,I,T,T,I);
cexp2=zeros(I,I,T,T,I);
%workw=zeros(I,I,T,T,I);
%workh=zeros(I,I,T,T,I);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
            [V,workh_, workw_, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2,V1,V2]=Vcouple(th,tw,jh,jw,i,PARREST,OUTC); %in period 1
             VC(jh,jw,th,tw,i)=V;
           % workh_check(jh,jw,th,tw,i)=workh; %out2(@() Vcouple(th,tw,jh,jw,i,p(i)));
           % workw_check(jh,jw,th,tw,i)=workw; %out3(@() Vcouple(th,tw,jh,jw,i,p(i)));
            
            if alloutput==1
                cexp1(jh,jw,th,tw,i)=conexp1; 
                cexp2(jh,jw,th,tw,i)=conexp2; 
                %workw(jh,jw,th,tw,i)=workw_; 
                %workh(jh,jw,th,tw,i)=workh_; 
            end
           
            Pw(jh,jw,th,tw,i,1,1)=Pnn; %out4(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,1,2)=Pnw0; %out5(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,1,3)=Pnw; %out6(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,2,1)=Pw0n; %out7(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,2,2)=Pw0w0; %out8(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,2,3)=Pw0w; %out9(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,3,1)=Pwn; %out10(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,3,2)=Pww0; %out11(@() Vcouple(th,tw,jh,jw,i,p(i)));
            Pw(jh,jw,th,tw,i,3,3)=Pww; %out12(@() Vcouple(th,tw,jh,jw,i,p(i)));
            VC1(jh,jw,th,tw,i)=V1; %out15(@() Vcouple(th,tw,jh,jw,i,p(i)));
            VC2(jh,jw,th,tw,i)=V2; %out16(@() Vcouple(th,tw,jh,jw,i,p(i))); % am Idefoing this unnecessarily?
                end
            end
        end
    end
end


% TODO: hear it should be possible to make 1 location twice as big (or
% something like that) exp(VC) in city times 2 should do it?
%{
PC=zeros(I,I,T,T,I);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                denom=sum(exp(VC(jh,jw,th,tw,:)/(sigmal)));
                for i=1:I
                PC(jh,jw,th,tw,i)=exp(VC(jh,jw,th,tw,i)/(sigmal))/denom;
                end
            end
        end
    end
end
%}
%{
cexp1(2,3,1,2,2)-cexp1(2,3,1,2,3) %typical offers - local offer is shit, so you would rather not.
cexp2(2,3,1,2,2)-cexp2(2,3,1,2,3) 

cexp1(3,2,1,2,2)-cexp1(3,2,1,2,3) %atypical offers - local offers is more enticing to take
cexp2(3,2,1,2,2)-cexp2(3,2,1,2,3) %it is surprising that this one is stronger effect for husbands? and it is not about lfp

% even without mm - matters more for wife to live close to her offer (not ic, offer),
% because they are more likely to take local offers regardless. so they
% pass on some nice matches if they are disadvantaged by commuting and insted take local offers. but this basic
% effect is tiny.

lambda=0.6;
lambda*(cexp1(2,3,1,2,2)-cexp1(2,3,1,2,3))+ (1-lambda)*(cexp2(2,3,1,2,2)-cexp2(2,3,1,2,3) )
exp(VC(2,3,1,2,:)./sigmal)/sum(exp(VC(2,3,1,2,:)./sigmal))
exp(VC(3,2,1,2,:)./sigmal)/sum(exp(VC(3,2,1,2,:)./sigmal)) % those that have atypical offers- more likely to prefer husband suburb
%}
doubleit=1;
if doubleit==1
PC=exp(VC./(repmat(sigmal*2,I,I,T,T,I)))./repmat(sum(exp(VC./(repmat(sigmal*2,I,I,T,T,I))),5),1,1,1,1,I);
else
    PC=exp(VC./(repmat(sigmal,I,I,T,T,I)))./repmat(sum(exp(VC./(repmat(sigmal,I,I,T,T,I))),5),1,1,1,1,I);
end

%{
DC1_=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:W
                        for ww=1:W
DC1_(jh,jw,th,tw,i,wh,ww)=   Pw(jh,jw,th,tw,i,wh,ww)*PC(jh,jw,th,tw,i)*NC(jh,jw,th,tw)/2;
                        end
                    end
                end
            end
        end
    end
end
%}
DC1=Pw.*repmat(PC,1,1,1,1,1,W,W).*repmat(NC,1,1,1,1,I,W,W)./repmat(2,I,I,T,T,I,W,W);

%{
ni=reshape(sum(sum(sum(sum(DC1,7),6))),T,T,I); % number of couples in th X tw X i
pjt=zeros(I,I,T,T); % conditional probability of jh jw given th and tw s
denoms=reshape(sum(sum(NC)),T,T);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                if denoms(th,tw)==0
                    pjt(jh,jw,th,tw)=0;
                else
                    pjt(jh,jw,th,tw)=NC(jh,jw,th,tw)/denoms(th,tw);
                end
            end
        end
    end
end
%}
pjt=NC./repmat(sum(sum(NC)),I,I);
pjt(isnan(pjt))=0;

%{
DC2=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I  
                    for wh=1:W
                        for ww=1:W
DC2(jh,jw,th,tw,i,wh,ww)=(mu^2*pjt(jh,jw,th,tw)*ni(th,tw,i) + ... % from both draw
                mu*(1-mu)*sum(pjt(jh,:,th,tw))*sum(sum(sum(DC1(:,jw,th,tw,i,:,:)))) + ... % just h draws
                (1-mu)*mu*sum(pjt(:,jw,th,tw))*sum(sum(sum(DC1(jh,:,th,tw,i,:,:)))) + ... % just w draws
                (1-mu)^2*sum(sum(DC1(jh,jw,th,tw,i,:,:))))*Pw(jh,jw,th,tw,i,wh,ww); % to noone draws - all stay put
            
b(jh,jw,th,tw,i,wh,ww)=sum(sum(sum(DC1(:,jw,th,tw,i,:,:)))); 
c(jh,jw,th,tw,i,wh,ww)=sum(sum(sum(DC1(jh,:,th,tw,i,:,:)))); 
                        end
                    end
                end
            end
        end
    end
end
%}

DC2=Pw.*((mu^2)*repmat(pjt,1,1,1,1,I,W,W).*repmat(sum(sum(sum(sum(DC1,7),6))),I,I,1,1,1,W,W)+...
    mu*(1-mu)*repmat(repmat(sum(pjt,2),1,I,1,1),1,1,1,1,I,W,W).*repmat(sum(sum(sum(DC1,6),7)),I,1,1,1,1,W,W)+...
    mu*(1-mu)*repmat(repmat(sum(pjt,1),I,1,1,1),1,1,1,1,I,W,W).*repmat(sum(sum(sum(DC1,6),7),2),1,I,1,1,1,W,W)+...
    ((1-mu)^2)*repmat(sum(sum(DC1,6),7),1,1,1,1,1,W,W)); %repmat(sum(sum(sum(DC1,6),7)),I,1,1,1,1,W,W) does not match exactly??
DC=DC1+DC2;


Hc=OUTC.('Hc') ; 
hC=zeros(I,I,T,T,I,W,W); 
%HCouple=zeros(I,I,T,T,I,W,W);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                for i=1:I
                    for wh=1:W
                        for ww=1:W
                            hi=wh>1;
                            wi=ww>1;
                            jh_=jh*(wh==3) + i*(wh==2|wh==1);
                            jw_=jw*(ww==3) + i*(ww==2|ww==1);
                            %lsh=hi*wi*lsb(jh_,jw_,i,1)+hi*(1-wi)*lsah(jh_,jw_,i); 
                            %lsw=hi*wi*lsb(jh_,jw_,i,2)+wi*(1-hi)*lsaw(jh_,jw_,i);
                            %Y=Ycc(th,tw,jh_,jw_,i,hi+1,wi+1);
                            hC(jh,jw,th,tw,i,wh,ww)=Hc(th,tw,jh_,jw_,i,3)*hi*wi+Hc(th,tw,jh_,jw_,i,1)*hi*(1-wi)+Hc(th,tw,jh_,jw_,i,2)*(1-hi)*wi; % hdC(Y,p(i)); % notice for some reason indeces have a different order

                        end 
                    end
                end
            end
        end
    end
end

HCouple=DC.*hC;


end


