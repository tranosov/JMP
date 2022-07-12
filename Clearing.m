%{
Couples and singles location

Function computing the difference between demand and supply in locations
%}



function clearing = Clearing(lp,EQS,PARREST)
%global D  AS  AC w uw uh us ch cw cs h hs alphaw betaw alphah betah mu lambda JLs Jw Jm HS Jc NC NS NSh NSw sigmaw 
%global WARNINGS %mm typeic HS NC NS
p=exp(lp);
%p(2)=1/p(2);
%p(3)=1/p(3);


HS=PARREST.('HS');


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

HSingleH=HSingle.*DSh ;
HSingleW=HSingle.*DSw ;
HSi=reshape(sum(sum(sum(HSingleH,4))),size(HS))+reshape(sum(sum(sum(HSingleW,4))),size(HS));
HCi=reshape(sum(sum(sum(sum(sum(sum(HCouple,7),6))))),size(HS)); 

clearing=(HSi+HCi-HS);

end