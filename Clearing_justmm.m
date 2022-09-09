%{
Couples and singles location

Function computing the difference between demand and supply in locations
%}



function [clearing] = Clearing_justmm(LA0,EQS,PARREST)
%global D  AS  AC uw uh us ch cw cs h hs alphaw betaw alphah betah mu lambda JLs Jw Jm HS Jc NC NS NSh NSw sigmaw 
global RESC %D Jw HS  %ssingle % WARNINGS 
if isempty(RESC)
    RESC=10^6;
end

HS=PARREST.('HS');
Jw=PARREST.('Jw');
sigmam=PARREST.('sigmam');
ssh_=PARREST.('sstaysingleh');
ssw_=PARREST.('sstaysinglew');
LA=LA0/(RESC);
wfh=PARREST.('wfh');

I=size(HS,2);
T=size(Jw,2);
if wfh>0
    T=T*2;
end
W=3;


% should I recalculate here the DC?
[~,~, DSh, DSw, VS,~]=DSingle(99,0,EQS,PARREST);
if (size(DSh,1)==1) 
    clearing=[1]*10^6; % penalty
    return
end
[DC, DC1, DC2,~,Pw,~,~,cexp1,cexp2,OUTC]=DCouple(99,1,LA,EQS,PARREST); 
if (size(DC,1)==1) 
    clearing=[1]*10^6; % penalty
    return
end

DSh_agr=sum(DSh,4);
DSw_agr=sum(DSw,4);
uhS=2*sum(sum(sum(DSh_agr.*VS)))./sum(sum(sum(DSh_agr)));
uwS=2*sum(sum(sum(DSw_agr.*VS)))./sum(sum(sum(DSw_agr)));

%global OUTC % when you call it here - you call the version created within the function, not before!
uh=OUTC.('uh') ;
uw=OUTC.('uw') ;
uhC_det=zeros(I,I,T,T,I,W,W); 
uwC_det=zeros(I,I,T,T,I,W,W);
uhC=zeros(I,I,T,T,I); 
uwC=zeros(I,I,T,T,I);
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
                           jh_=jh*(wh==3) + i*(wh==2|wh==1);
                           jw_=jw*(ww==3) + i*(ww==2|ww==1);
                           if wwi>0
                           uhC_det(jh,jw,th,tw,i,wh,ww)=uh(th,tw,jh_,jw_,i,wwi) ; 
                           uwC_det(jh,jw,th,tw,i,wh,ww)=uw(th,tw,jh_,jw_,i,wwi) ;
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

DC_agr=sum(sum(DC,7),6);
DC1_agr=sum(sum(DC1,7),6);
DC2_agr=sum(sum(DC2,7),6);

uhC=sum(sum(sum(sum(sum(DC1_agr.*uhC)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uhC)))))./sum(sum(sum(sum(sum(DC2_agr)))));
uwC=sum(sum(sum(sum(sum(DC1_agr.*uwC)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uwC)))))./sum(sum(sum(sum(sum(DC2_agr)))));
uhC=round(uhC,6);
uwC=round(uwC,6);
uhS=round(uhS,6);
uwS=round(uwS,6);



% no - need to do it separately for DC1 and DC2! check if same?

% Notice - the distribution does not depend on the share that gets married
% here! (given lambda of course)

%sigmam=1;
M=sum(sum(sum(DSh_agr)))*(1/(1/3 + ssh_*2/3)) ; %+ sum(sum(sum(sum(sum(DC_agr)))));
F=sum(sum(sum(DSw_agr)))*(1/(1/3 + ssw_*2/3)) ; %+ sum(sum(sum(sum(sum(DC_agr)))));

clmm=  F*(exp((uwC-uwS)/sigmam)/(1+exp((uwC-uwS)/sigmam)))  - M*(exp((uhC-uhS)/sigmam)/(1+exp((uhC-uhS)/sigmam)));

%ssh=1/(1+exp((uhC-uhS)/sigmam)) ;
%ssw=1/(1+exp((uwC-uwS)/sigmam)) ;
clearing=[clmm*(1/3)];

%todo - make a version where lambda is being identified separately - and
%instead of here changing lambda - change

end