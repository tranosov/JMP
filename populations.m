function [DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC, DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS,EQS,PARREST] = populations(p,LA,EQS,PARREST,FORCEFIT)


HS=PARREST.('HS');
Jw=PARREST.('Jw');
sigmam=PARREST.('sigmam');
ssh_=PARREST.('sstaysingleh');
ssw_=PARREST.('sstaysinglew');

I=size(HS,2);
T=size(Jw,2);
W=3;


[DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC]=DCouple(p,1,LA,EQS,PARREST);
[DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS]=DSingle(p,1,EQS,PARREST);


DSh_agr=sum(DSh,4);
DSw_agr=sum(DSw,4);
uhS=2*sum(sum(sum(DSh_agr.*VS)))./sum(sum(sum(DSh_agr)));
uwS=2*sum(sum(sum(DSw_agr.*VS)))./sum(sum(sum(DSw_agr)));

%global OUTC % when you call it here - you call the version created within the function, not before!
uh=OUTC.('uh') ;
uw=OUTC.('uw') ;
uhC_det=zeros(I,I,T,T,I,W,W); 
uwC_det=zeros(I,I,T,T,I,W,W);
uhC_=zeros(I,I,T,T,I); 
uwC_=zeros(I,I,T,T,I);
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
uhC_(jh,jw,th,tw,i)= sum(sum(Pw(jh,jw,th,tw,i,:,:).*uhC_det(jh,jw,th,tw,i,:,:))) + cexp1(jh,jw,th,tw,i);
uwC_(jh,jw,th,tw,i)= sum(sum(Pw(jh,jw,th,tw,i,:,:).*uwC_det(jh,jw,th,tw,i,:,:))) + cexp2(jh,jw,th,tw,i);
                
                end
            end
        end
    end
end

DC_agr=sum(sum(DC,7),6);
DC1_agr=sum(sum(DC1,7),6);
DC2_agr=sum(sum(DC2,7),6);

uhC=sum(sum(sum(sum(sum(DC1_agr.*uhC_)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uhC_)))))./sum(sum(sum(sum(sum(DC2_agr)))));
uwC=sum(sum(sum(sum(sum(DC1_agr.*uwC_)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uwC_)))))./sum(sum(sum(sum(sum(DC2_agr)))));
% no - need to do it separately for DC1 and DC2! check if same?

% Notice - the distribution does not depend on the share that gets married
% here! (given lambda of course)

M=sum(sum(sum(DSh_agr)))*(1/(1/3 + ssh_*2/3)) ; %+ sum(sum(sum(sum(sum(DC_agr)))));
F=sum(sum(sum(DSw_agr)))*(1/(1/3 + ssw_*2/3)) ; %+ sum(sum(sum(sum(sum(DC_agr)))));

if FORCEFIT==1 % recalibrate THETA and THETAHW (for within estimation routine)
    params=PARREST.('params');
    THETA=params{'THETA',:};  
    THETAHW=params{'THETAHW',:};  
    dTHETAHW=-(uhC-uhS) - log( (M/F)*(1+exp(-(uwC-uwS)/sigmam))  -1)*sigmam; % making marriage market clear
    params{'THETAHW',:}=THETAHW + dTHETAHW;
    % make women match
    if exist('momentall','var')
        smarried=momentall{'smarried',:};
    else
        smarried=1-ssw_;
    end
    dTHETA= log((smarried)/(1-smarried))*sigmam - (uwC-uwS);
    params{'THETA',:}=THETA + dTHETA;
    
    PARREST.('params')=params;
    uhC=uhC+dTHETA+dTHETAHW;
    uwC=uwC+ dTHETA ;  % does this work?
    OUTC.('uh')=uh+(dTHETA+dTHETAHW)/2 ; % per period!
    OUTC.('uw')= uw+ (dTHETA)/2  ;
       
end


if FORCEFIT==0.5 % only THETAHW, not THETA - just so I can treat lambda as a parameter.
    params=PARREST.('params');
    %THETA=params{'THETA',:};  
    THETAHW=params{'THETAHW',:};  
    dTHETAHW=-(uhC-uhS) - log( (M/F)*(1+exp(-(uwC-uwS)/sigmam))  -1)*sigmam; % making marriage market clear
    params{'THETAHW',:}=THETAHW + dTHETAHW;
    PARREST.('params')=params;
    uhC=uhC+dTHETA+dTHETAHW;
    uwC=uwC+ dTHETA ;  % does this work?
    OUTC.('uh')=uh+(dTHETAHW)/2 ; % per period!
    OUTC.('uw')= uw  ;
       
end

clmm=  F*(exp((uwC-uwS)/sigmam)/(1+exp((uwC-uwS)/sigmam)))  - M*(exp((uhC-uhS)/sigmam)/(1+exp((uhC-uhS)/sigmam)));
if clmm>10^(0)
    fprintf('Marrriage market uncleared - on purpose?')
    FORCEFIT=FORCEFIT
    clmm=clmm
end
ssh=1/(1+exp((uhC-uhS)/sigmam)) ;
ssw=1/(1+exp((uwC-uwS)/sigmam)) ;
%log((1-0.15)/0.15)*sigmam = dU

CONSTSH=(ssh*M*(2/3)+(1/3)*M  )/(sum(sum(sum(DSh_agr))));
CONSTSW=(ssw*F*(2/3)+(1/3)*F  )/(sum(sum(sum(DSw_agr))));
CONSTC=(1/2)*(M*(2/3)*(1-ssh)+F*(2/3)*(1-ssw))/(sum(sum(sum(sum(sum(DC_agr)))))) ;
% it should be true that M*(2)*(1-ssh) = F*(2)*(1-ssw)

HSingleH=HSingle.*DSh*CONSTSH ;
HSingleW=HSingle.*DSw*CONSTSW ;

HSingle=HSingleH+HSingleW;
HCouple=(HCouple)*CONSTC; % todo: make a function for all these sums

DSh=DSh*CONSTSH; % has to be gender specific!
DSw=DSw*CONSTSW;
DS=DSh+DSw;
DC=DC*CONSTC;
DC1=DC1*CONSTC;
DC2=DC2*CONSTC;

% issue: marriage market does not clear here, if it is not forced to. I
% should have resolved.

end