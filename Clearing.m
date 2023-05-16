%{
Couples and singles location

Function computing the difference between demand and supply in locations
%}



function clearing = Clearing(lp,EQS,PARREST)
%global D  AS  AC w uw uh us ch cw cs h hs alphaw betaw alphah betah mu lambda JLs Jw Jm HS Jc NC NS NSh NSw sigmaw 
%global WARNINGS %mm typeic HS NC NS
global show
p=exp(lp);
%p(2)=1/p(2);
%p(3)=1/p(3);


HS=PARREST.('HS');
ssh_=PARREST.('sstaysingleh');
ssw_=PARREST.('sstaysinglew');

[DS,HSingle,DSh, DSw]=DSingle(p,1,EQS,PARREST);
if (size(DS,1)==1)
    clearing=[1,1,1]*10^6; % penalty
    return
end

[DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC]=DCouple(p,1,99,EQS,PARREST);
if  (size(DC,1)==1) 
    clearing=[1,1,1]*10^6; % penalty
    return
end
% (DSi+DCi-HS); IF COUPLES OCCUPY JUST 1 UNIT OF HOUSING!
%measure=(clearing.^2)./HS;

    DSh_agr=sum(DSh,4);
    DSw_agr=sum(DSw,4);
    DC_agr=sum(sum(DC,7),6);
    
    M=sum(sum(sum(DSh_agr)))+(sum(sum(sum(sum(sum(DC_agr)))))); %sum(sum(sum(DSh_agr)))*(1/(1/3 + ssh_*2/3)) ; %+ sum(sum(sum(sum(sum(DC_agr)))));
    F=sum(sum(sum(DSw_agr)))+(sum(sum(sum(sum(sum(DC_agr)))))); %sum(sum(sum(DSw_agr)))*(1/(1/3 + ssw_*2/3)) ; %+ sum(sum(sum(sum(sum(DC_agr)))));
    CONSTSH=(ssh_*M*(2/3)+(1/3)*M  )/(sum(sum(sum(DSh_agr))));
    CONSTSW=(ssw_*F*(2/3)+(1/3)*F  )/(sum(sum(sum(DSw_agr))));
    
    CONSTC=(1/2)*(M*(2/3)*(1-ssh_)+F*(2/3)*(1-ssw_))/(sum(sum(sum(sum(sum(DC_agr)))))) ;
    % it should be true that M*(2)*(1-ssh) = F*(2)*(1-ssw)

HSingleH=HSingle.*DSh*CONSTSH ;
HSingleW=HSingle.*DSw*CONSTSW ;
HSi=reshape(sum(sum(sum(HSingleH,4))),size(HS))+reshape(sum(sum(sum(HSingleW,4))),size(HS));
HCi=reshape(sum(sum(sum(sum(sum(sum(HCouple,7),6))))),size(HS))*CONSTC ; % todo: make a function for all these sums

%HSingleH=HSingle.*DSh ;
%HSingleW=HSingle.*DSw ;
%HSi=reshape(sum(sum(sum(HSingleH,4))),size(HS))+reshape(sum(sum(sum(HSingleW,4))),size(HS));
%HCi=reshape(sum(sum(sum(sum(sum(sum(HCouple,7),6))))),size(HS)); 

clearing=(HSi+HCi-HS);

end