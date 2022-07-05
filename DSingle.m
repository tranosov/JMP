%{
Couples and singles location

Demand functions for singles.

Job offer in j X type t.
%}

function [DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS]=DSingle(p,alloutput,EQS,PARREST)
global WARNINGS 

D=PARREST.('D');
sigmal=PARREST.('sigmal');
NS=PARREST.('NS');
NSh=PARREST.('NSh');
NSw=PARREST.('NSw');


% Create housind demand equations in each location:
I=size(D,1);
T=size(NS,2);
W=3;
WARNINGS=0; %reset
OUTS=laborsS(p,alloutput,EQS,PARREST);

if WARNINGS>0
    DS=999*10^6;
    HSingle=999;
    DSh=999;
    DSw=999;
    VS=999;
    Pws=999;
    cexps=999;
    OUTS=999;
    return
end

%Vsingle(1,2,1,1)

%Vsingle(1,1,2,1) % commute to city from my suburb
%Vsingle(1,3,2,1) % commute to worse opposite suburb? - set as a little bit
%Vsingle(1,1,3,1) % commute to city from bad suburb - should be MOST commutes
%Vsingle(1,2,1,1) % commute to my suburb from city (both good) - many!
%Vsingle(1,3,1,1) % commute to a bad suburb from a city



Pws=zeros(I,T,I,W);
VS=zeros(I,T,I); %values
cexps=zeros(I,T,I);
for j=1:I
    for t=1:T
        for i=1:I
            [VS_,works,Pn,Pw0,Pw_,conexp]=Vsingle(t,j,i,PARREST,OUTS);
            Pws(j,t,i,1)=Pn; %out3(@() Vsingle(t,j,i,p(i)));
            Pws(j,t,i,2)=Pw0; %out4(@() Vsingle(t,j,i,p(i)));
            Pws(j,t,i,3)=Pw_; %out5(@() Vsingle(t,j,i,p(i)));
            VS(j,t,i)=VS_; %Vsingle(t,j,i,p(i));
            cexps(j,t,i)=conexp;
           
        end
    end
end

%{
PS=zeros(I,T,I); %probabilities
for j=1:I
    for t=1:T
        denom=sum(exp(VS(j,t,:)/sigmal));
        for i=1:I
            PS(j,t,i)=exp(VS(j,t,i)/sigmal)/denom;
        end
    end
end
%}
PS=exp(VS./repmat(sigmal,I,T,I))./repmat(sum(exp(VS./repmat(sigmal,I,T,I)),3),1,1,I);

%{
DS=zeros(I,T,I,W); % singles in categories
for j=1:I
    for t=1:T
        for i=1:I
            for w=1:W
                DS(j,t,i,w)=PS(j,t,i)*NS(j,t)*Pw(j,t,i,w);
            end
        end
    end
end
%}
DS=Pws.*repmat(NS,1,1,I,W).*repmat(PS,1,1,1,W);
DSh=Pws.*repmat(NSh,1,1,I,W).*repmat(PS,1,1,1,W);
DSw=Pws.*repmat(NSw,1,1,I,W).*repmat(PS,1,1,1,W);

%{
DSh=zeros(I,T,I,W); 
for j=1:I
    for t=1:T
        for i=1:I
            for w=1:W
                DSh(j,t,i,w)=PS(j,t,i)*NSh(j,t)*Pw(j,t,i,w);
            end
        end
    end
end

DSw=zeros(I,T,I,3);
for j=1:I
    for t=1:T
        for i=1:I
            for w=1:W
                DSw(j,t,i,w)=PS(j,t,i)*NSw(j,t)*Pw(j,t,i,w);
            end
        end
    end
end
%}

Hs=OUTS.('Hs') ;

%HSingle=zeros(I,T,I,w); 
hS=zeros(I,T,I,W); 
%commuteS=zeros(I,T,I,W);
worksS=zeros(I,T,I,W);
for j=1:I
    for t=1:T
        for i=1:I
            for w=1:W
                worksS(j,t,i,w)=(w>1);
                %Y=(w==3)*Yss(t,j,i) + (w==2)*Yss(t,i,i);
                %(w==3)*lss(j,i)*(w>1)*w1(lss(j,i))+ (w==2)*lss(i,i)*(w>1)*w1(lss(i,i));
                hS(j,t,i,w)=(w==3)*Hs(t,j,i) + (w==2)*Hs(t,i,i); %hdS(Y,p(i));
                %commuteS(j,t,i,w)=D(i,j)*(w==3) + D(i,i)*(w==2);
            end
        end
    end
end
HSingle=hS; %DS.*

%commuteS_raw=sum(sumS(commuteS.*DS.*worksS))./sum(sumS(DS.*worksS));

% restructure?
% j offered and j used would still not be enough, because having i=j used
% can be for 2 different reasons which would matter for welfar. 
% could do j offered and 0 1 2 for 'working'. This is most 'model', but not
% at all 'data', but still probably makes sense.
% Still I guess it is ok to aggregate them a bit? 
%PwS=Pw;
end

%todo: single men vs women!
