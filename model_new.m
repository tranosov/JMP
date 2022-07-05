
% additions: piw_(dh,dw), H/2

function [EQS,PARREST]=...
    model_new(toinputs,AS1,A1,AC_sub, AC2,dmin,dminS,d, a,ah,aw,b,bh,bw,timeh,timew,l,Jcenter, Jschool,...
    Wcenter, Wschool,Hcenter, Hschool,Jseg,sw,sel,wa,wb, wlinear, cneces,...
    ces,ceh,cew,leffecth,leffectw,eta,etaC,muprob,plocal_,sstaysingleh,sstaysinglew,mA_,mI_,kappa_,...
    gamd,PHID,deltaw,PHI,rhod,PHIW_,line,muw,muwL,swL,crra,crrat,crrah, crrax,pih0_,piw0_,pi_,pish_,pisw_,piel_,pitheta_,qel_,...
    d21,d31,Sup1,Sup2,Sup3,SupS,addS,out,daddS,aaddS,pophous,...
    XI,wc, PID,typeic_,pid_,wgap_raw,mm_,lsbar, NLY,THETA,THETAHW,sigmam,MtoF) 
rng(357)
global VERBOSE
mm=mm_;
typeic=typeic_;
%todo: redo as a structure!

%todo: create pih_,piw_ as functions of d? figure out how to make the
%dependence not obviously gendered but turn out to be gendered?

if wlinear==1 % derivative of w(ls)*ls
    w1=@(ls,ic) wa+wb.*ls+ic.*wc;
    w2=w1; %wage*wgap; % notice - for now wb and wa are gender neutral!
    w1_d=@(ls,ic)wa+2*wb.*ls+ic.*wc; % derivative of ls*w(ls)=w(ls)+ls*w'(ls)
    w2_d=@(ls,ic)wa+2*wb.*ls+ic.*wc;
elseif wlinear==2 % exp
    w1=@(ls,ic) wa.*exp(ls*wb+ic*wc);
    w2=@(ls,ic) wa.*exp(ls*wb+ic*wc+wgap_raw); %wage*wgap; % notice - for now wb and wa are gender neutral!
    w1_d=@(ls,ic)wa.*exp(ls*wb+ic*wc).*(1+ls*wb);
    w2_d=@(ls,ic)wa.*exp(ls*wb+ic*wc+wgap_raw).*(1+ls*wb);
end
subsub=sqrt(d21-1+line^2/4)+sqrt(d31-1+line^2/4);
D=[0,1*d21,1*d31;...
    1*d21,0,subsub;...
    1*d31,subsub,0]; % line between 1-triangle and 2-line
Dinsub=[0,0,0; 0,1,0; 0,0,1];
D=D.*d+dmin.*(D==0)+(dminS-dmin).*Dinsub;
AS=[A1+AS1,0,0]; % amenities as percieved by singles. city better.
AC=[A1+0,AC_sub+AC2,AC_sub]; % city vs suburb is no different, but one suburb is a better school district
if addS==1
    if out==1 % the new neighborhood is put on the outside to be 'as far as possible' from both other suburbs (triangle)
        newsub=sqrt(line^2/4 + (daddS + sqrt(1-line^2/4))^2);
        D=[0,1,1,daddS;... #lets assume d21=d31=1!
            1,0,1*line,newsub;...
            1,1*line,0,newsub;
            daddS,  newsub,newsub,0];
        Dinsub=[0,0,0,0; 0,1,0,0; 0,0,1,0;0,0,0,0]; % treat is as a suburb? no!
        D=D.*d+dmin.*(D==0)+(dminS-dmin).*Dinsub;
        if aaddS==1 % something in between
            AS=[AS,0.5*AS(1)+0.5*AS(3)];
            AC=[AC,0.5*AC(1)+0.5*AC(3)];
        elseif aaddS==2 % suburb
            AS=[AS,AS(3)];
            AC=[AC,AC(3)];
        elseif aaddS==3 % city
            AS=[AS,AS(1)];
            AC=[AC,AC(1)];
        elseif aaddS==4 % sliding!
            AS=[AS,newsub*AS(1)+(1-newsub)*AS(3)];
            AC=[AC,newsub*AC(1)+(1-newsub)*AC(3)];
        elseif aaddS==23 % proportion as original
            AS=[AS,(1/3)*AS(1)+(2/3)*AS(3)];
            AC=[AC,(1/3)*AC(1)+(2/3)*AC(3)];
        else
            AS=[AS,0];
            AC=[AC,0];
        end
    end
end

%eta=1;

% make share spend on housing independent
alphas=a;
betas=b;
alphah=a*ah;
betah=b*bh;
alphaw=a*aw;
betaw=b*bw;
tih=timeh;
tiw=timew;
sigmal=sel; % use in creating proababilities only!

%+ const*(1-betas*d) could make sense. pi will have to be much more gendered. Issue - for d to matter, const will be potentially high, which could make lambda matter VERY DIFFERENTLY


mu=muprob; % probability of switching jobs.
plocal=plocal_;
lambda=l; %bargaining weight - HUSBANDS BARGAINING WEIGHT?

% Distributions of location of job offers by type! matrix iXt. Columns sum
% up to 1. First column is the more 'female' industry.
jl1=Jcenter/(2+Jcenter);
jl3=Jschool*(1-jl1)/(1+Jschool);
jl2=1-jl1-jl3;
JL11=Wcenter*jl1;
JL12=Hcenter*jl1;
JL31=(jl3/jl2)*Wschool*(1-JL11)/(1+(jl3/jl2)*Wschool);
JL32=(jl3/jl2)*Hschool*(1-JL12)/(1+(jl3/jl2)*Hschool);
JL21=1-JL11-JL31;
JL22=1-JL12-JL32;


JLs=[ JL11, JL12 ;... % center
      JL21, JL22;... % suburb
      JL31, JL32];
      %ones(3,2)*1/3; 
      %;[ 1,0.4 ;... % center
      %0, 0.3;... % suburb
      %0, 0.3]; %suburb with a good school district
% This will have to be matched in a complicated way, because I observe only the distribution of jobs taken.
if addS==1
    Jcenter_=2*(JL11+JL12)/(JL21+JL22+JL31+JL32);
    JLs=[JLs*(Jcenter_+2)/(Jcenter_+3);
        1/(Jcenter_+3),1/(Jcenter_+3)];
        
end

T=size(JLs);
T=T(2);
% Distribution of men/women in sectors. SECOND sector is more womanly.
Jw=[1-Jseg,Jseg]; %[0.3,0.7]; % check in the data. Later swich to more groups (adding a vertical dimension)
Jm=[Jseg,1-Jseg];
%Jc=Jw'*Jm;
Z_=[1,1;1,1];
Jc=logitmatching(Jw,Jm,Z_); % rows men, columns women


% Supply of housing
Sup2_=Sup2+SupS/2; % actually this should be one or the other!
Sup3_=Sup3+SupS/2;
N=size(D);
I=N(1);
HS=[1000*(1+Sup1/100),1000*(1+Sup2_/100),1000*(1+Sup3_/100)]; 
if I==4
    HS=[1000,1000,1000,1000];
end
% 1500 singles and 1500 couples: 750 couples + 750 couples=1500+1500*2=4500 people
% So there are ultimately twice as many people in couples.

% Implies number of single people in jxt (job offer location X industry)
if pophous==1
    NN=sum(HS);
else
    NN=I*1000;
end
ssingleh=1/3+sstaysingleh*2/3; % CRAZY NUMBER???? 0.43 -INT THAT TOO MUCH?
ssinglew=1/3+sstaysinglew*2/3;
sm=1/(1+1/MtoF);
NSh=NN*ssingleh*sm*JLs.*repmat(Jm,I,1); %sm~0.5
NSw=NN*ssinglew*(1-sm)*JLs.*repmat(Jw,I,1);
NS=NSh+NSw; %NSh/NSw should be M/F here

% Implied number of couples in jxjxtxt
NC=zeros(I,I,T,T);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
            NC(jh,jw,th,tw)=Jc(th,tw)*JLs(jh,th)*JLs(jw,tw); %shares
            end
        end
    end
end

NC=NN*(1-(ssingleh+ssinglew)/2)*NC*0.5; %m&w together (so half as many couples as there are allocated housing units to couples
%NC(jh,jw,th,tw)

mIL=muwL-sqrt(3)*swL;
mAL=muwL+sqrt(3)*swL;
mI=muw-sqrt(3)*sw;
mA=muw+sqrt(3)*sw;



%pih_=@(dh,dw) pih0_; %@(d) pih0_*((1-betas*d)^PID);
piw_=@(dh,dw)  min(1,max(piw0_ + betah*(dh-dw)*pid_,0)); %min(0.999,max(piw0_ + betah*(dh-dw)*pid_,0.0001)); %=@(d) piw0_*((1-betas*d)^PID);

if crra==1
    uu=@(x) max(log(x),-1.0000e+10);
else
    uu=@(x) max(((x).^(1-crra)-1)./(1-crra),-1.0000e+10); % the higher crra the more concave
end
if crrat==1
    uutime=@(x) max(log(x),-1.0000e+10);
else
    uutime=@(x) max(((x).^(1-crrat)-1)./(1-crrat),-1.0000e+10); % the higher crra the more concave
end
if crrax==1
    uux=@(x) max(log(x),-1.0000e+10);
else
    uux=@(x) max(((x).^(1-crrax)-1)./(1-crrax),-1.0000e+10); % the higher crra the more concave
end
if crrah==1
    uuh=@(x) max(log(x),-1.0000e+10); % notice - rival but non-excludable!
else
    uuh=@(x) max(((x).^(1-crrah)-1)./(1-crrah),-1.0000e+10); % the higher crra the more concave
end

%xw=@(ls,d) ((piw_(d))^(1/crra)/((piw_(d))^(1/crra)+(1-lambda)^(1/crra)))*(1-ls*alphaw - d*betaw+tiw); % it is missing lew, if I wanted a variation there
%xh=@(ls,d) ((pih_(d))^(1/crra)/((pih_(d))^(1/crra)+(lambda)^(1/crra)))*(1-ls*alphah - d*betah+tih); % it is missing leh, if I wanted a variation there
%((piw_(d))^(1/crra)/((piw_(d))^(1/crra)+(1-lambda)^(1/crra))) - percent of
%your time spent on housework?

PGT_=@(xh,xw,dh,dw) (piw_(dh,dw)==1)*((xw.^(1-piel_)).^(1/(1-piel_))) + ...
    (piw_(dh,dw)==0)*((xh.^(1-piel_)).^(1/(1-piel_))) + ...
    (piw_(dh,dw)>0)*(piw_(dh,dw)<1)*(max(piw_(dh,dw),0.000000001).*xw.^(1-piel_)+ (max(1-piw_(dh,dw),0.000000001)).*xh.^(1-piel_) ).^(1/(1-piel_)); % pi_.*(piw_.*xw.^(1-piel_)+ (1-piw_).*xh.^(1-piel_) ).^(1/(1-piel_)); 
PG_= @(xh,xw,q,dh,dw) (PGT_(xh,xw,dh,dw)^(1-qel_)+(q*pitheta_)^(1-qel_))^(1/(1-qel_));
PG=@(xh,xw,q,dh,dw) pi_.*uux(PG_(xh,xw,q,dh,dw));
PGsh=@(xh,q) pish_.*uux((xh^(1-qel_)+(q*pitheta_)^(1-qel_))^(1/(1-qel_))); %pish_.*xh;
PGsw=@(xw,q)  pisw_.*uux((xw^(1-qel_)+(q*pitheta_)^(1-qel_))^(1/(1-qel_)));

%{
fdistw=@(dw)PHIW_*betas*dw;
if rhod==1
    fdist=@(dh,dw)((betah*dh)^(1-deltaw))*((dw*betaw)^deltaw); %(1/(deltaw^deltaw + (1-deltaw)^(1-deltaw)))
elseif rhod==0
    fdist=@(dh,dw)max(deltaw,1-deltaw)*min(betah*dh/(1-deltaw),betaw*dw/deltaw);     
else
    fdist=@(dh,dw)((((1-deltaw)^(1/rhod))/((1-deltaw)^(1/rhod)+deltaw^(1/rhod)))*(betah*dh)^((rhod-1)/rhod)+...
        ((deltaw^(1/rhod))/((1-deltaw)^(1/rhod)+deltaw^(1/rhod)))*(dw*betaw)^((rhod-1)/rhod))^(rhod/(rhod-1)); 
end
fdistS=@(d)betas*d; % to keep things in the same scale - make it in time. And avoid 0
%}

if gamd<1000
    Fdist=@(dh,dw)((((1-deltaw)^(gamd))/((1-deltaw)^(gamd)+deltaw^(gamd)))*(1-betah*dh)^(gamd)+...
        ((deltaw^(gamd))/((1-deltaw)^(gamd)+deltaw^(gamd)))*(1-dw*betaw)^(gamd))^(1/gamd); 
else
    Fdist=@(dh,dw) 2*max( (((1-deltaw))/((1-deltaw)+deltaw))*(1-betah*dh), ((deltaw)/((1-deltaw)+deltaw))*(1-dw*betaw) ); 
end
% infty~ max?
FdistS=@(d)(1-betas*d);

uw=@(work,dh,dw,a,c,h,ls,xh_,xw_,q,ic) (cew*uu(-cneces+c)*(c >0)+...
    (-1.0000e+10)*(-cneces+c <=0) + a+ etaC*uuh(h/2) ...
    +PHID*(Fdist(dh,dw)-1)...
    +PG(xh_,xw_,q,dh,dw) ...
    +leffectw*uutime((1-work*ls*alphaw - work*dw*betaw+tiw-xw_))... 
    +XI*ls*ic)+ THETA/2;
    %-PHI*fdist(dh,dw)+...
    %-fdistw(dw)+...
uh=@(work,dh,dw,a,c,h,ls,xh_,xw_,q,ic) (ceh*uu(-cneces+c)*(c >0)+...
    (-1.0000e+10)*(-cneces+c <=0) + a+ etaC*uuh(h/2)...
    +PHID*(Fdist(dh,dw)-1) ...
    +PG(xh_,xw_,q,dh,dw) ...
    +leffecth*uutime((1-work*ls*alphah - work*dh*betah-tih-xh_))...  
    +XI*ls*ic)+ THETA/2 + THETAHW/2;
    %-PHI*fdist(dh,dw)+...
    %-fdistw(dw
V=@(workh,workw,dh,dw,a,ch,cw,h,lsh,lsw,xh,xw,q,ich,icw) lambda*uh(workh,dh,dw,a,ch,h,lsh,xh,xw,q,ich)+(1-lambda)*uw(workw,dh,dw,a,cw,h,lsw,xh,xw,q,icw);   
V_lambda=@(workh,workw,dh,dw,a,ch,cw,h,lsh,lsw,xh,xw,q,ich,icw,lambda) lambda*uh(workh,dh,dw,a,ch,h,lsh,xh,xw,q,ich)+(1-lambda)*uw(workw,dh,dw,a,cw,h,lsw,xh,xw,q,icw);   

us=@(work,d,a,c,h,ls,x,q,ic)(ces*uu(-cneces+c)*(c >0)+...
    (-1.0000e+10)*(-cneces+c <=0) + a+ eta*uuh(h)+ leffecth*uutime((1-work*ls*alphas - x-work*d*betas -tih))...
    +PHID*(FdistS(d)-1)...
    +PGsh(x,q) ...
    +XI*ls*ic);
    %-PHI*fdistS(d) ...
ush=@(work,d,a,c,h,ls,x,q,ic)(ces*uu(cneces+c)*(c >0)+...
    (-1.0000e+10)*(cneces+c <=0) + a+ eta*uuh(h)+ leffecth*uutime((1-work*ls*alphas - x-work*d*betas -tih))...
    +PHID*(FdistS(d)-1)...
    +PGsh(x,q) ...
    +XI*ls*ic);
usw=@(work,d,a,c,h,ls,x,q,ic)(ces*uu(-cneces+c)*(c >0)+...
    (-1.0000e+10)*(-cneces+c <=0) + a+ eta*uuh(h)+ leffecth*uutime((1-work*ls*alphas - x- work*d*betas -tih))...
    +PHID*(FdistS(d)-1)...
    +PGsw(x,q) ...
    +XI*ls*ic);

%XX=((1-lambda)/lambda)^(1/crra);
%Z=(lambda^(1/crra)+(1-lambda)^(1/crra))^crra;
Ys=@(ls,ic) ls.*w1(ls,ic) + 0.5*NLY ; %workh*w+workw*w*wgap;
Yc=@(lsh, lsw,ich,icw) lsh.*w1(lsh,ich)+...
    lsw.*w2(lsw,icw) + NLY; % using wage functions

 
% other variables - as functions of prices and multiplier!
hdS=@(mu,p) (eta/(p*mu))^(1/crrah); %(Y)*(p+(p*ceh/eta)^(1/crra))^(-1); % assumes ceh=cew
hdC=@(mu,p)  ( 2^(1-1/crrah))*(etaC/(p*mu))^(1/crrah); %(Y)*(p+(lambda^(1/crra)+(1-lambda)^(1/crra))*(p*ceh/etaC)^(1/crra))^(-1);

ch=@(mu,lambda)  (lambda*ceh/mu)^(1/crra)+cneces; %Y*lambda^(1/crra)./((lambda^(1/crra)+(1-lambda)^(1/crra)) + ((etaC/ceh)^(1/crra))*...
cw=@(mu,lambda)  ((1-lambda)*cew/mu)^(1/crra)+cneces; %Y*(1-lambda)^(1/crra)./((lambda^(1/crra)+(1-lambda)^(1/crra)) + ((etaC/ceh)^(1/crra))*...
cs=@(mu)  (ces/(mu))^(1/crra) +cneces;  %$Y*((1 + (eta/ceh)^(1/crra)*p^(1-1/crra))^(-1));
cc=@(mu,lambda)  ((lambda)^(1/crra)+(1-lambda)^(1/crra) )*(ceh/mu)^(1/crra) +2*cneces;

if pitheta_>0
    qs=@(mu,x) 0; %(mu/(pitheta_*pish_*x.^((1-crrax)*(1-pitheta_)))).^(1/(pitheta_*(1-crrax)-1)); % would have to change if sh != sw
    qc=@(mu,xh,xw) 0; %(mu/(pitheta_*pi_*PGT_(xh,xw).^((1-crrax)*(1-pitheta_)))).^(1/(pitheta_*(1-crrax)-1));
else
    qs=@(mu,x) 0;
    qc=@(mu,xh,xw) 0;
end


CONS=1000;
% prepare equations to solve for mu (only in spacial cases it has a closed form)
multS_eq=@(Y,p,mu,x) [(- Y + cs(mu) + hdS(mu,p) + qs(mu,x))]; % will solve a system of equations for (mu,x)
multC_eq=@(Y,p,mu,xh,xw,lambda) CONS*[(- Y + cc(mu,lambda) + hdC(mu,p) + qc(mu,xh,xw))];

%multS=@(ls,p,ic) (1/ceh)*(Ys(ls,ic).^crra).*(1 + p^(1-1/crra)*(eta/ceh)^(1/crra) )^(-crra); %1/mult!
%multC=@(lsh,lsw,p,ich,icw) (1/ceh)*((Yc(lsh, lsw,ich,icw)).^crra).*(Z^(1/crra) + p^(1-1/crra)*(etaC/ceh)^(1/crra) )^(-crra);
%I think it is ok that the 1/mult for couples is twice as big (their
%multiplier is smaller, because each person is 'a half'. The leisure also counts as half.

% singles 
Leissh= @(mu,x) (1./((1-pitheta_)*pish_*(qs(mu,x).^(pitheta_*(1-crrax)))*x.^((1-pitheta_)*(1-crrax)-1) )).^(1/crrat); %note - this is not leisure this is leisure-tih
mulSh=@(mu,x,ic) alphas.*leffecth./(Leissh(mu,x).^crrat) -XI*ic; % - der wrt h
lssh=@(mu,x,d) 1-betah*d -x - Leissh(mu,x)-tih;
lssh_eq=  @(mu,x,p,d,ic) (mulSh(mu,x,ic) - mu*w1_d(lssh(mu,x,d) ,ic))*(lssh(mu,x,d)<lsbar)+(lssh(mu,x,d)==1)*0 + (lssh(mu,x,d)>lsbar)*(CONS)*(lssh(mu,x,d)-lsbar)^1;
%mulSw=@(ls,d,ic) alphas.*leffecth./((1-alphas.*ls -betas*d+tih -xsw_fun(ls,d)).^crrat) -XI*ic; % - der wrt h
%lssw_eq=  @(ls,p,d,ic) multS(ls,p,ic).*mulSw(ls,d,ic) - w1_d(ls,ic); % to solve for ls (d)


% couples - both work
%Leish_= @(q,xh,xw) (    lambda./( (1-pitheta_)*pi_*(           q^(pitheta_*(1-crrax)))*(PGT_(xh,xw)^((1-pitheta_)*(1-crrax)-1))*(pih_*(xh^(-piel_))*(PGT_(xh,xw)^piel_)) ) ).^(1/crrat);  

Leish= @(mu,xh,xw,dh,dw,lambda) (    lambda./( pi_*PGT_(xh,xw,dh,dw)^(-crrax)*(max(1-piw_(dh,dw),0.0000001))*...
    (xh^(-piel_))*(PGT_(xh,xw,dh,dw)^piel_) ) ).^(1/crrat);  
%(    lambda./( (1-pitheta_)*pi_*(qc(mu,xh,xw)^(pitheta_*(1-crrax)))*(PGT_(xh,xw,dh,dw)^((1-pitheta_)*(1-crrax)-1))*((max(1-piw_(dh,dw),0.0000001))*...
%    (xh^(-piel_))*(PGT_(xh,xw,dh,dw)^piel_)) ) ).^(1/crrat);  
Leisw= @(mu,xh,xw,dh,dw,lambda) (    (1-lambda)./( pi_*PGT_(xh,xw,dh,dw)^(-crrax)*(max(piw_(dh,dw),0.0000001))*...
    (xw^(-piel_))*(PGT_(xh,xw,dh,dw)^piel_) ) ).^(1/crrat);  
%((1-lambda)./( (1-pitheta_)*pi_*(qc(mu,xh,xw)^(pitheta_*(1-crrax)))*(PGT_(xh,xw,dh,dw)^((1-pitheta_)*(1-crrax)-1))*(max(piw_(dh,dw),0.0000001)*...
%    (xw^(-piel_))*(PGT_(xh,xw,dh,dw)^piel_)) ) ).^(1/crrat); 
mulh=@(mu,xh,xw,dh,dw,ic,lambda) lambda*alphah.*leffecth./((Leish(mu,xh,xw,dh,dw,lambda)).^crrat) -XI*ic;
mulw=@(mu,xh,xw,dh,dw,ic,lambda) (1-lambda)*alphaw.*leffectw./((Leisw(mu,xh,xw,dh,dw,lambda)).^crrat) -XI*ic;

lsh=@(mu,o2,o3,dh,dw,lambda) (piw_(dh,dw)>0 && piw_(dh,dw)<1 )*(max(0,1-betah*dh -o2 - Leish(mu,o2,o3,dh,dw,lambda)-tih)) + (piw_(dh,dw)==1)*(1-betah*dh -o2)+ (piw_(dh,dw)==0)*(1-betah*dh -o2 -Leish(mu,o2,0.000001,dh,dw,lambda)-tih);
lsw=@(mu,o2,o3,dh,dw,lambda) (piw_(dh,dw)>0 && piw_(dh,dw)<1 )*(max(0,1-betaw*dw -o3 - Leisw(mu,o2,o3,dh,dw,lambda)-tiw)) + (piw_(dh,dw)==0)*(1-betah*dh -o3)+ (piw_(dh,dw)==1)*(1-betah*dh -o3 -Leisw(mu,0.000001,o3,dh,dw,lambda)-tiw); 

lshY=@(mu,o2,o3,dh,dw,lambda) (piw_(dh,dw)>0 && piw_(dh,dw)<1 )*(1-betah*dh -o2 - Leish(mu,o2,o3,dh,dw,lambda)-tih) + (piw_(dh,dw)==1)*(1-betah*dh -o2)+ (piw_(dh,dw)==0)*(1-betah*dh -o2 -Leish(mu,o2,0.000001,dh,dw,lambda)-tih);
lswY=@(mu,o2,o3,dh,dw,lambda) (piw_(dh,dw)>0 && piw_(dh,dw)<1 )*(1-betaw*dw -o3 - Leisw(mu,o2,o3,dh,dw,lambda)-tiw) + (piw_(dh,dw)==0)*(1-betah*dh -o3)+ (piw_(dh,dw)==1)*(1-betah*dh -o3 -Leisw(mu,0.000001,o3,dh,dw,lambda)-tiw); 


xh_fun=@(mu,o2,o3,dh,dw) (piw_(dh,dw)>0 && piw_(dh,dw)<1 )*o2 + (piw_(dh,dw)==1)*0 + (piw_(dh,dw)==0)*o2;
xw_fun=@(mu,o2,o3,dh,dw) (piw_(dh,dw)>0 && piw_(dh,dw)<1 )*o3 + (piw_(dh,dw)==1)*o3 + (piw_(dh,dw)==0)*0;


lsb_eq1= @(mu,xh,xw,p,dh,dw,ich,icw,lambda) (mulh(mu,xh,xw,dh,dw,ich,lambda) - w1_d(lsh(mu,xh,xw,dh,dw,lambda),ich)*mu)*(lsh(mu,xh,xw,dh,dw,lambda)<lsbar) +...
    (lsh(mu,xh,xw,dh,dw,lambda)>lsbar)*(CONS)*(lsh(mu,xh,xw,dh,dw,lambda)-lsbar)^1 +(lsh(mu,xh,xw,dh,dw,lambda)<=0)*10^6;
lsb_eq2= @(mu,xh,xw,p,dh,dw,ich,icw,lambda) (mulw(mu,xh,xw,dh,dw,icw,lambda) - w2_d(lsw(mu,xh,xw,dh,dw,lambda),ich)*mu)*(lsw(mu,xh,xw,dh,dw,lambda)<lsbar) + ...
(lsw(mu,xh,xw,dh,dw,lambda)>lsbar)*(CONS)*(lsw(mu,xh,xw,dh,dw,lambda)-lsbar)^1 + (lsw(mu,xh,xw,dh,dw,lambda)<=0)*10^6;
lsb_eq = @(mu,xh,xw,p,dh,dw,ich,icw,lambda)[ lsb_eq1(mu,xh,xw,p,dh,dw,ich,icw,lambda),lsb_eq2(mu,xh,xw,p,dh,dw,ich,icw,lambda)];

% couples - one works    
% just husband works
lsah_eq= @(mu,xh,xw,p,dh,dw,ic,lambda) [(mulh(mu,xh,xw,dh,dw,ic,lambda)  - w1_d(lsh(mu,xh,xw,dh,dw,lambda),ic)*mu)*(lsh(mu,xh,xw,dh,dw,lambda)<lsbar) + (lsh(mu,xh,xw,dh,dw,lambda)>lsbar)*(CONS)*(lsh(mu,xh,xw,dh,dw,lambda)-lsbar)^1, ...
    (Leisw(mu,xh,xw,dh,dw,lambda)  - (1-xw)+tiw)*CONS];
%xx=xh/xw

% just wife working
lsaw_eq= @(mu,xh,xw,p,dh,dw,ic,lambda) [(mulw(mu,xh,xw,dh,dw,ic,lambda)  - w2_d(lsw(mu,xh,xw,dh,dw,lambda),ic)*mu)*(lsw(mu,xh,xw,dh,dw,lambda)<lsbar) + (lsw(mu,xh,xw,dh,dw,lambda)>lsbar)*(CONS)* (lsw(mu,xh,xw,dh,dw,lambda)-lsbar)^1 , ...
    (Leish(mu,xh,xw,dh,dw,lambda)  - (1-xh)+tih)*CONS];
%xx=xw/xh

% for piw=0 or piw=1: 
lsb_eq_xh0= @(mu,Lh,xw,p,dh,dw,ich,icw,lambda) [(lambda*alphah.*leffecth./((Lh-tih).^crrat) - XI*ich - w1_d(1-betah*dh - Lh,ich)*mu)*(1-betah*dh - Lh<lsbar)  + (1-betah*dh - Lh <=0)*10^6 + (1-betah*dh - Lh >lsbar)*(CONS)*(1-betah*dh - Lh -lsbar)^1,...
    (mulw(mu,0,xw,dh,dw,icw,lambda) - w2_d(lsw(mu,0,xw,dh,dw,lambda),ich)*mu)*(lsw(mu,0,xw,dh,dw,lambda)<lsbar) + (lsw(mu,0,xw,dh,dw,lambda)<=0)*10^6 + (lsw(mu,0,xw,dh,dw,lambda)>lsbar)*(CONS)*(lsw(mu,0,xw,dh,dw,lambda)-lsbar)^1];
lsb_eq_xw0= @(mu,xh,Lw,p,dh,dw,ich,icw,lambda) [(mulh(mu,xh,0,dh,dw,ich,lambda) - w1_d(lsh(mu,xh,0,dh,dw,lambda),ich)*mu)*(lsh(mu,xh,0,dh,dw,lambda)<lsbar) +  (lsh(mu,xh,0,dh,dw,lambda)<=0)*10^6 + (lsh(mu,xh,0,dh,dw,lambda)>lsbar)*(CONS)*(lsh(mu,xh,0,dh,dw,lambda)-lsbar)^1, ...
    ((1-lambda)*alphah.*leffecth./((Lw-tiw).^crrat) - XI*icw - w2_d(1-betah*dw - Lw,icw)*mu)*(1-betah*dw - Lw<lsbar) + (1-betah*dw - Lw<=0)*10^6 + (1-betah*dw - Lw>lsbar)*(CONS)*(1-betah*dw - Lw - lsbar)^1];

lsah_eq_xh0= @(mu,Lh,xw,p,dh,dw,ic,lambda) [(lambda*alphah.*leffecth./((Lh-tih).^crrat)  - XI*ic - w1_d(1-betah*dh - Lh,ic)*mu)*(1-betah*dh - Lh<lsbar) + (1-betah*dh - Lh>lsbar)*(CONS)*(1-betah*dh - Lh-lsbar)^1 + (1-betah*dh - Lh<=0)*10^6 ,...
    (Leisw(mu,0,xw,dh,dw,lambda)  - (1-xw)+tiw)*CONS];
lsah_eq_xw0= @(mu,xh,xw,p,dh,dw,ic,lambda) [(mulh(mu,xh,0,dh,dw,ic,lambda)  - w1_d(lsh(mu,xh,0,dh,dw,lambda),ic)*mu)*(lsh(mu,xh,0,dh,dw,lambda)<lsbar) + (lsh(mu,xh,0,dh,dw,lambda)>lsbar)*(CONS)*(lsh(mu,xh,0,dh,dw,lambda) - lsbar)^1 + (lsh(mu,xh,0,dh,dw,lambda)<=0)*10^6,...
    xw*CONS];

lsaw_eq_xh0= @(mu,xh,xw,p,dh,dw,ic,lambda) [(mulw(mu,0,xw,dh,dw,ic,lambda)  - w2_d(lsw(mu,0,xw,dh,dw,lambda),ic)*mu)*(lsw(mu,0,xw,dh,dw,lambda)<lsbar) + (lsw(mu,0,xw,dh,dw,lambda)>lsbar)*(CONS)*(lsw(mu,0,xw,dh,dw,lambda) - lsbar)^1 + (lsw(mu,0,xw,dh,dw,lambda)<=0)*10^6, ...
    xh*CONS];
lsaw_eq_xw0= @(mu,xh,Lw,p,dh,dw,ic,lambda) [ ((1-lambda)*alphah.*leffecth./((Lw-tiw).^crrat) -XI*ic - w2_d(1-betah*dw - Lw,ic)*mu)*(1-betah*dw - Lw<lsbar)+ (1-betah*dw - Lw>lsbar)*(CONS)*(1-betah*dw - Lw-lsbar)^1 + (1-betah*dw - Lw<=0)*10^6, ...
    (Leish(mu,xh,0,dh,dw,lambda)  - (1-xh)+tih)*CONS];



if toinputs==1
    options = optimoptions('fsolve','MaxIter',500,'MaxFunctionEvaluations',500,...
            'FunctionTolerance',10^(-6),'Display','off','Algorithm','levenberg-marquardt',...
            'StepTolerance', 10^(-10),'ScaleProblem','jacobian'); %'trust-region'); %'StepTolerance',10^(-15),'ScaleProblem','jacobian',
    x0=0.05;
    mu0=ces*(1/( Ys(0.23,0)-1))^(crra); %0.7402; % figure out!
    p=1;
    inputs0S=zeros(T,I,I,2);

    for j=1:I
        for i=1:I
            for t=1:T
                d=D(i,j);
                [~,~,~,~,~,low]=matchdist(i,j,t,0,0,0,0,typeic_,D,mm,JLs,betah);
                ic=1-low;
                fn=@(mu,x) [lssh_eq(mu,x,p,d,ic) ,multS_eq(Ys(lssh(mu,x,d),ic),p,mu,x)];
                fn=@(in) fn(in(1),in(2));
                [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,[mu0,x0],options); 
                output_=output;
                if (~isreal(output)) | (output(1)<=0) | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                        if VERBOSE
                            fprintf('INIT: Warning on singles')
                        end
                        [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,[mu0,x0].*sign(fn([mu0,x0])).*[-1,1]*1.1,options); 
                         if (~isreal(output)) | (output(1)<=0) | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                            if VERBOSE
                                fprintf('INIT: Warning on singles - 2')
                            end
                            
                            [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,[mu0,x0].*sign(fn([mu0,x0])).*[-1,1]*1.0005,options); 
                             if (~isreal(output)) | (output(1)<=0) | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                fprintf('INIT: Warning on singles - 3')
                                output=output_;
                             else
                                %fprintf('Solved');   
                             end
                         else
                             %fprintf('Solved');    
                         end
                end
                inputs0S(t,j,i,:)=real(output);

            end
        end
    end
    ich0=0.0154;
    icw0=0.0154;
    mu0=  (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(0.26,0.2,ich0,icw0)-2.1))^(crra); %l*(1/(w1_d(0.258,0)))*(1/(0.69)^crrat); %0.42; %0.3902;0.6856
    mu00= (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(0.2672,0,ich0,icw0)-1.4))^(crra); %l*(1/(w1_d(0.26,0)))*(1/(0.6856)^crrat) ; %0.45; %0.3885; 
    mu000=  (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(0,0.2672,ich0,icw0)-1.4))^(crra); %l*(1/(w1_d(0.26,0)))*(1/(0.6856)^crrat) ; %0.45; %0.3885; 
    
    inputs0=zeros(T,T,I,I,I,3,3);
    for i=1:I
        for jh=1:I
            for jw=1:I 
                for th=1:T
                    for tw=1:T
                         dh=D(i,jh);
                         dw=D(i,jw);
                         inputs0(th,tw,jh,jw,i,1,:)=[mu00,0.03*100,0.1890*100];
                         inputs0(th,tw,jh,jw,i,2,:)=[mu000,0.13*100,0.06*100];
                         inputs0(th,tw,jh,jw,i,3,:)=[mu0,0.04*100,0.11*100];
                    end
                end
            end
        end
    end

    inputs=zeros(T,T,I,I,I,3,3);

    p0=[1.1,1,1];    

    for jh=1:I
        for jw=1:I 
            for th=1:T
                for tw=1:T
                    for i=1:I
                        p=p0(i);
                        dh=D(i,jh);
                        dw=D(i,jw);
                        [~,~,~,~,~,low]=matchdist(i,jh,th,0,0,0,0,typeic,D,mm,JLs,betah);
                        ich=1-low;
                        [~,~,~,~,~,low]=matchdist(i,jw,tw,0,0,0,0,typeic,D,mm,JLs,betah);
                        icw=1-low;
                  
                        if piw_(dh,0)>0 && piw_(dh,0)<1 
                            fn=@(mu,xh,xw) [lsah_eq(mu,xh,xw,p,dh,0,ich,lambda) ,multC_eq(Yc(lsh(mu,xh,xw,dh,0,lambda),0,ich,icw),p,mu,xh,xw,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output1,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,1,:),options) ;
                        end
                        if piw_(dh,0)==1
                            %fprintf('xh=0');
                            fn=@(mu,Lh,xw) [lsah_eq_xh0(mu,Lh,xw,p,dh,0,ich,lambda) ,multC_eq(Yc(1-betah*dh - Lh,0,ich,icw),p,mu,0,xw,lambda)];
                            fn=@(x)fn(x(1),x(2)/100,x(3)/100);
                            [output1,~,EXITFLAG]=fsolve(fn,[inputs0(th,tw,jh,jw,i,1,1),(1-betah*dh-0.26)*100,inputs0(th,tw,jh,jw,i,1,3)],options) ;
                        end
                        if piw_(dh,0)==0
                            %fprintf('xw=0');
                            fn=@(mu,xh,xw) [lsah_eq_xw0(mu,xh,xw,p,dh,0,ich,lambda) ,multC_eq(Yc(lsh(mu,xh,xw,dh,0,lambda),0,ich,icw),p,mu,xh,xw,lambda)]; % one equation is just xw=0
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output1,~,EXITFLAG]=fsolve(fn,[inputs0(th,tw,jh,jw,i,1,1),inputs0(th,tw,jh,jw,i,1,2),0],options) ;
                        end
                        output1_=output1;
                        if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))

                            %fprintf('Warning');
                            %cl=fn(output2);
                            inp_=(piw_(0,dw)>0 && piw_(0,dw)<1)*reshape(inputs0(th,tw,jh,jw,i,1,:),[1,3])*0.9 + ...
                                (piw_(0,dw)==1)*[inputs0(th,tw,jh,jw,i,1,1),(1-betah*dh-0.26)*100,inputs0(th,tw,jh,jw,i,1,3)]*0.9+...
                                (piw_(0,dw)==0)*[inputs0(th,tw,jh,jw,i,1,1),inputs0(th,tw,jh,jw,i,1,2),(1-betaw*dw)*100]*0.9;
                            
                            [output1,~,EXITFLAG]=fsolve(fn,inp_,options); % to lazy, do better
                            
                            if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                if VERBOSE
                                    fprintf('Warning 2');
                                end
                                inp_=(piw_(0,dw)>0 && piw_(0,dw)<1)*reshape(inputs0(th,tw,jh,jw,i,1,:),[1,3])*1.1 + ...
                                (piw_(0,dw)==1)*[inputs0(th,tw,jh,jw,i,1,1),(1-betah*dh-0.26)*100,inputs0(th,tw,jh,jw,i,1,3)]*1.1+...
                                (piw_(0,dw)==0)*[inputs0(th,tw,jh,jw,i,1,1),inputs0(th,tw,jh,jw,i,1,2),(1-betaw*dw)*100]*1.1;
                                [output1,~,EXITFLAG]=fsolve(fn,inp_,options); % to lazy, do better
                                if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                    fprintf('INIT: Warning w0')
                                    output1=output1_;
                                else
                                    %fprintf('Solved');
                                end
                            else
                                %fprintf('Solved');
                            end
                        end
                        
                        inputs(th,tw,jh,jw,i,1,:)=real(output1);

                        
                        if piw_(0,dw)>0 && piw_(0,dw)<1 
                            fn=@(mu,xh,xw) [lsaw_eq(mu,xh,xw,p,0,dw,icw,lambda) ,multC_eq(Yc(0,lsw(mu,xh,xw,0,dw,lambda),ich,icw),p,mu,xh,xw,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output2,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,2,:),options); 
                        end
                        if piw_(0,dw)==1
                            %fprintf('xh=0');
                            fn=@(mu,xh,xw) [lsaw_eq_xh0(mu,xh,xw,p,0,dw,icw,lambda) ,multC_eq(Yc(0,lsw(mu,xh,xw,0,dw,lambda),ich,icw),p,mu,xh,xw,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output2,~,EXITFLAG]=fsolve(fn,[inputs0(th,tw,jh,jw,i,2,1),0,inputs0(th,tw,jh,jw,i,2,3)],options); 
                        end
                        if piw_(0,dw)==0
                            %fprintf('xw=0');
                            fn=@(mu,xh,Lw) [lsaw_eq_xw0(mu,xh,Lw,p,0,dw,icw,lambda) ,multC_eq(Yc(0,1-betah*dw - Lw,ich,icw),p,mu,xh,0,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output2,~,EXITFLAG]=fsolve(fn,[inputs0(th,tw,jh,jw,i,2,1),inputs0(th,tw,jh,jw,i,2,2),1-betaw*dw-0.26],options); 
                        end
                        output2_=output2;
                        if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                            if VERBOSE
                            fprintf('Warning');
                            end
                            %cl=fn(output2);
                            inp_=(piw_(0,dw)>0 && piw_(0,dw)<1)*reshape(inputs0(th,tw,jh,jw,i,2,:),[1,3])*0.9 + ...
                                (piw_(0,dw)==1)*[inputs0(th,tw,jh,jw,i,2,1),1-betah*dh-0.26,inputs0(th,tw,jh,jw,i,2,3)]*0.9+...
                                (piw_(0,dw)==0)*[inputs0(th,tw,jh,jw,i,2,1),inputs0(th,tw,jh,jw,i,2,2),1-betaw*dw-0.23]*0.9;
                            
                            [output2,~,EXITFLAG]=fsolve(fn,inp_,options); % to lazy, do better
                            
                            if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                fprintf('INIT: Warning 0w')
                                output2=output2_;
                            else
                                %fprintf('Solved');
                            end
                            
                        end
                        inputs(th,tw,jh,jw,i,2,:)=real(output2); 


                        if piw_(dh,dw)>0 && piw_(dh,dw)<1
                            fn=@(mu,xh,xw) [lsb_eq(mu,xh,xw,p,dh,dw,ich,icw,lambda),...
                                multC_eq(Yc(lsh(mu,xh,xw,dh,dw,lambda),lsw(mu,xh,xw,dh,dw,lambda),ich,icw),p,mu,xh,xw,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output3,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,3,:),options);
                        end
                        if piw_(dh,dw)==1
                            fprintf('xh=0');
                            fn=@(mu,Lh,xw) [lsb_eq_xh0(mu,Lh,xw,p,dh,dw,ich,icw,lambda),...
                                multC_eq(Yc(1-Lh-betah*dh,lsw(mu,0,xw,dh,dw,lambda),ich,icw),p,mu,0,xw,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                           [output3,~,EXITFLAG]=fsolve(fn,[inputs0(th,tw,jh,jw,i,3,1),(1-betah*dh-0.26)*100,inputs0(th,tw,jh,jw,i,3,3)],options);
                        end
                        if piw_(dh,dw)==0
                            fprintf('xw=0');
                            fn=@(mu,xh,Lw) [lsb_eq_xw0(mu,xh,Lw,p,dh,dw,ich,icw,lambda),...
                                multC_eq(Yc(lsh(mu,xh,0,dh,dw,lambda),1-betah*dw - Lw,ich,icw),p,mu,xh,0,lambda)];
                            fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                            [output3,~,EXITFLAG]=fsolve(fn,[inputs0(th,tw,jh,jw,i,3,1),inputs0(th,tw,jh,jw,i,3,2),(1-betaw*dw-0.23)*100],options);
                        end
                        output3_=output3;
                        if ~isreal(output3) | output3(1)<=0 | output3(2)<=0 |  lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0  |...
                                lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                            if VERBOSE
                                fprintf('Warning');
                            end
                            %options = optimoptions('fsolve','MaxIter',5000,'MaxFunctionEvaluations',5000,...
                            %            'FunctionTolerance',10^(-8),'Display',iter_,...
                            %            'StepTolerance', 10^(-12),'Algorithm','trust-region');
                            const=0.85;
                            inp_=(piw_(dh,dw)>0 && piw_(dh,dw)<1)*reshape(inputs0(th,tw,jh,jw,i,3,:),[1,3])*const + ...
                                (piw_(dh,dw)==1)*[inputs0(th,tw,jh,jw,i,3,1),(1-betah*dh-0.26)*100,inputs0(th,tw,jh,jw,i,3,3)]*0.85+...
                                (piw_(dh,dw)==0)*[inputs0(th,tw,jh,jw,i,3,1),inputs0(th,tw,jh,jw,i,3,2),(1-betaw*dw-0.23)*100]*0.85;
                            [output3,~,EXITFLAG]=fsolve(fn, inp_,options); 
                            %todo: if really no solution - somehow code
                            %this option never to be used?
                            if ~isreal(output3) | output3(1)<=0 | output3(2)<=0 | lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0  |...
                                lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                if VERBOSE
                                    fprintf('Warning2');
                                end
                                const=0.7;
                                inp_=(piw_(dh,dw)>0 && piw_(dh,dw)<1)*reshape(inputs0(th,tw,jh,jw,i,3,:),[1,3])*const + ...
                                    (piw_(dh,dw)==1)*[inputs0(th,tw,jh,jw,i,3,1),1-betah*dh-0.26,inputs0(th,tw,jh,jw,i,3,3)]*0.85+...
                                    (piw_(dh,dw)==0)*[inputs0(th,tw,jh,jw,i,3,1),inputs0(th,tw,jh,jw,i,3,2),1-betaw*dw-0.23]*0.85;
                                [output3,~,EXITFLAG]=fsolve(fn, inp_,options); 
                                if ~isreal(output3) | output3(1)<=0 | output3(2)<=0 | lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0  |...
                                    lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                    fprintf('INIT: Warning ww')
                                    output3=output3_;
                                else
                                    %fprintf('Solved');
                                end
                            else
                                %fprintf('Solved');
                            end
     
                        end 
                      inputs(th,tw,jh,jw,i,3,:)=real(output3);
              
                    end
                end
            end
        end
    end

end




EQS = struct();

EQS.('lssh_eq') = lssh_eq;   
EQS.('multS_eq') = multS_eq;    
EQS.('multC_eq') = multC_eq; 
EQS.('lsah_eq') = lsah_eq; 
EQS.('lsaw_eq') = lsaw_eq; 
EQS.('lsb_eq') = lsb_eq; 
EQS.('lsah_eq_xh0') = lsah_eq_xh0; 
EQS.('lsaw_eq_xh0') = lsaw_eq_xh0; 
EQS.('lsb_eq_xh0') = lsb_eq_xh0; 
EQS.('lsah_eq_xw0') = lsah_eq_xw0; 
EQS.('lsaw_eq_xw0') = lsaw_eq_xw0; 
EQS.('lsb_eq_xw0') = lsb_eq_xw0; 

EQS.('lssh') = lssh; 
EQS.('Ys') = Ys; 
EQS.('hdS') =  hdS;
EQS.('cs') = cs;
EQS.('Yc') = Yc; 
EQS.('hdC') =  hdC;
EQS.('ch') = ch;
EQS.('cw') = cw;
EQS.('cc') = cc;
EQS.('us') = us;
EQS.('lsh') = lsh;
EQS.('lsw') = lsw;

EQS.('V_nolambda') = V;
EQS.('V') = V_lambda; % can change the weighting
EQS.('uh') = uh;
EQS.('uw') = uw;
EQS.('xh_fun') = xh_fun;
EQS.('xw_fun') = xw_fun;
EQS.('utilcs') = @(x)ces*uu(x) ;
EQS.('utills') = @(ls,x,d) uutime(1-x-ls*alphas - (ls>0)*d*betas-tih) ;
EQS.('utilxs') =  PGsh;
EQS.('utilcls') = @(d)PHID*(FdistS(d)-1);
EQS.('utilhs') = @(x)eta*uuh(x) ;


EQS.('utilc') = @(x)ceh*uu(x) ;
EQS.('utill') = @(ls,x,d) uutime(1-x-ls*alphas - (ls>0)*d*betas-tih) ;
%EQS.('utilcw') = @(x)cew*uu(x) ;
%EQS.('utillw') = @(ls,x,d) uutime(1-x-ls*alphas - (ls>0)*d*betas+tih) ;
EQS.('utilho') = @(x)etaC*uuh(x/2) ;
EQS.('utilx') =  PG;
EQS.('utilcl') = @(dh,dw) PHID*(Fdist(dh,dw)-1);
EQS.('leis') = @(ls,x,d) 1-x-ls*alphas - (ls>0)*d*betas-tih ;

%todo: think -  singles vs couples 'public goods': do I implicitely have
%couples value it twice as much? (I think so) this does not mean the demand
%is twice as big btw. since everybody has 'more' of it, they can spend
%less!

EQS.('w1') = w1;
EQS.('w2') = w2;
EQS.('piw') = piw_;
EQS.('uXi')= @(ls,ic) XI*ls*ic;

EQS.('mu0') = @(ls,l) 0.5*(1/(w1_d(ls,0)))*(1/(l)^crrat); %ignore
if toinputs==1
    EQS.('inputsS') = real(inputs0S);
    EQS.('inputs') = real(inputs);
else
    EQS.('inputsS') = 0;
    EQS.('inputs') = 0;
end

PARREST.('lambda')=lambda;
PARREST.('ceh')=ceh;
PARREST.('ces')=ces;
PARREST.('crra')=crra;
PARREST.('JLs')=JLs;
PARREST.('mm')=mm;
PARREST.('plocal')=plocal;
PARREST.('mu')=mu;
PARREST.('betah')=betah;
PARREST.('mA')=mA;
PARREST.('mI')=mI;
PARREST.('mAL')=mAL;
PARREST.('mIL')=mIL;
PARREST.('NC')=NC;
PARREST.('NS')=NS;
PARREST.('NSh')=NSh;
PARREST.('NSw')=NSw;
PARREST.('AS')=AS;
PARREST.('AC')=AC;
PARREST.('Jw')=Jw;
PARREST.('HS')=HS;
PARREST.('sigmam')=sigmam;
PARREST.('sigmal')=sigmal;
PARREST.('D')=D;
PARREST.('typeic')=typeic;
PARREST.('XI')=XI;
PARREST.('Jc')=Jc;
PARREST.('sstaysingleh')=sstaysingleh;
PARREST.('sstaysinglew')=sstaysinglew;

end



 
