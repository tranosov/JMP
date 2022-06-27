%{
Couples and singles location

Function computing the value for singles.

t - type
j - job offer location
i - residence location
p - price of housing

Maybe later add wage, now treat as a global variable!
%}

function [V,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2,V1,V2]=Vcouple(th,tw,jh,jw,i, PARREST,OUTC)

JLs=PARREST.('JLs');
mu=PARREST.('mu');
D=PARREST.('D');

[V1,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2]=Vcouplesharesmatch_2(th,tw,jh,jw,i, PARREST,OUTC);

% prepare all options of jh/jw (3x3 matrix)
N=size(D);
EV2=zeros(N);
for ii=1:N(1)
    for jj=1:N(2)
        EV2(ii,jj)=Vcouplesharesmatch_2(th,tw,ii,jj,i, PARREST,OUTC);
    end
end

EV2w=EV2(jh,:)*JLs(:,tw);
EV2h=JLs(:,th)'*EV2(:,jw);
EV2both=JLs(:,th)'*EV2*JLs(:,tw);
V2=mu^2*EV2both + mu*(1-mu)*EV2h + (1-mu)*mu*EV2w + (1-mu)*(1-mu)*Vcouplesharesmatch_2(th,tw,jh,jw,i, PARREST,OUTC); % here I am recomputing stuff unnecessarily - I should have created all of them somewhere and only use them


% so far work/consumption decision in period 1 do not affecrt V2. In that
% case the total value is just a sum of two independent values.
% And the first value is the same function as the second value. 

V=V1+V2;
end


