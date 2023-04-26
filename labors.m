function OUTC=labors(pp,alloutput,lambda,EQS,PARREST)
%global EQS JLs AC typeic D betah WARNINGS mm params
global WARNINGS VERBOSE
global IN


OUTC = struct();

typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
AC=PARREST.('AC');
betah=PARREST.('betah');
ceh=PARREST.('ceh');
crra=PARREST.('crra');
wfh_share=PARREST.('wfh_share');

N=size(JLs);
I=N(1);
T=N(2);
if pp==99
    params=PARREST.('params');
    pp=[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ]; 
    if I==4
    pp=[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:},params{'p0_4',:} ]; 
    end
end

if lambda==99
    params=PARREST.('params');
    lambda=params{'LA0',:}; 
end

lsah_eq=EQS.('lsah_eq');  
lsaw_eq=EQS.('lsaw_eq');  
lsb_eq=EQS.('lsb_eq'); 
lsah_eq_xh0=EQS.('lsah_eq_xh0');  
lsaw_eq_xh0=EQS.('lsaw_eq_xh0');  
lsb_eq_xh0=EQS.('lsb_eq_xh0'); 
lsah_eq_xw0=EQS.('lsah_eq_xw0');  
lsaw_eq_xw0=EQS.('lsaw_eq_xw0');  
lsb_eq_xw0=EQS.('lsb_eq_xw0'); 

multC_eq=EQS.('multC_eq');
lsh=EQS.('lsh') ;
lsw=EQS.('lsw');
ch_fun=EQS.('ch') ;
cw_fun=EQS.('cw');
Yc=EQS.('Yc');
utilh=EQS.('uh');
utilw=EQS.('uw');
utilc=EQS.('utilc');
utill=EQS.('utill');
%utilcw=EQS.('utilcw');
utilx=EQS.('utilx');
utilho=EQS.('utilho');
utilcl=EQS.('utilcl');
hdC=EQS.('hdC');
V=EQS.('V');
piw_=EQS.('piw');
uXi=EQS.('uXi');
mu0_=EQS.('mu0');
xh_fun=EQS.('xh_fun');
xw_fun=EQS.('xw_fun');

wfh=PARREST.('wfh');
T0=T;
if wfh>0
    T=T0*2;
end
hs_wfh=PARREST.('hs_wfh');


if alloutput
    lh=zeros(T,T,I,I,I,3);
    lw=zeros(T,T,I,I,I,3);
    xh=zeros(T,T,I,I,I,3);
    xw=zeros(T,T,I,I,I,3);
    mus=zeros(T,T,I,I,I,3);
    Y=zeros(T,T,I,I,I,3);
    ch=zeros(T,T,I,I,I,3);
    cw=zeros(T,T,I,I,I,3);
    uh=zeros(T,T,I,I,I,3);
    uw=zeros(T,T,I,I,I,3);
    uch=zeros(T,T,I,I,I,3);
    ucw=zeros(T,T,I,I,I,3);
    ulh=zeros(T,T,I,I,I,3);
    ulw=zeros(T,T,I,I,I,3);
    ux=zeros(T,T,I,I,I,3);
    uho=zeros(T,T,I,I,I,3);
    ucl=zeros(T,T,I,I,I,3);
    piw_dd=zeros(T,T,I,I,I,3);
    uXih=zeros(T,T,I,I,I,3);
    uXiw=zeros(T,T,I,I,I,3);
    ich_=zeros(T,T,I,I,I,3);
    icw_=zeros(T,T,I,I,I,3);
    %leih=zeros(T,T,I,I,I,3);
    %leiw=zeros(T,T,I,I,I,3);
end
%inputs=zeros(T,T,I,I,I,3,3);
vc=zeros(T,T,I,I,I,3); % last - only h, only w, both work
Hc=zeros(T,T,I,I,I,3);
TOL=10^(-15);
STEPTOL=10^(-10); % would need more iterations! plus by melo byt <TOL!

% lingering issue: still I have sometimes different guesses reading to
% slightly differenc cl. And persistently sow.
% I think this is being stuck at the boundary of lsbar!

% but more troubelingly - sometimes I get the sense that this does not have
% a unique solution!

iter_='off'; %'iter'
options = optimoptions('fsolve','MaxIter',2000,'MaxFunctionEvaluations',2000,...
                 'FunctionTolerance',TOL,'Display',iter_,'Algorithm','trust-region','StepTolerance', STEPTOL); %better for when your guess is already good
             
options0 = optimoptions('fsolve','MaxIter',1000,'MaxFunctionEvaluations',1000,...
        'FunctionTolerance',TOL,'Display',iter_,'Algorithm','trust-region',...
        'StepTolerance', STEPTOL); %,'StepTolerance', STEPTOL); %'trust-region') ;%,...
        %'StepTolerance', 10^(-20)); %'ScaleProblem','jacobian'
        %options = optimoptions('fsolve','MaxIter',5000,'MaxFunctionEvaluations',5000,...
         %   'FunctionTolerance',TOL,'Display','off','Algorithm','levenberg-marquardt'); 
         % it seems to me that 'levenberg-marquardt' is just faster?
             
options1 = optimoptions('fsolve','MaxIter',5000,'MaxFunctionEvaluations',5000,...
                 'FunctionTolerance',TOL,'Display',iter_,'Algorithm','levenberg-marquardt','StepTolerance', STEPTOL);         
         
ich0=0.0135; % also adjust! 
icw0=0.0135;
xh0_2=0.03;
xh0_3=0.13;
xh0_1=0.04;
xw0_2=0.189;
xw0_3=0.06;
xw0_1=0.11;

mu00= (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(lsh(1,xh0_2,xw0_2,10,8,lambda),0,ich0,icw0)-1.4))^(crra); %l*(1/(w1_d(0.26,0)))*(1/(0.6856)^crrat) ; %0.45; %0.3885; 
mu000=  (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(0,lsw(1,xh0_3,xw0_3,10,8,lambda),ich0,icw0)-1.4))^(crra); %l*(1/(w1_d(0.26,0)))*(1/(0.6856)^crrat) ; %0.45; %0.3885; 
mu0=  (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(lsh(1,xh0_1,xw0_1,10,8,lambda),lsw(1,xh0_1,xw0_1,10,8,lambda),ich0,icw0)-2.1))^(crra); %l*(1/(w1_d(0.258,0)))*(1/(0.69)^crrat); %0.42; %0.3902;0.6856

for ii=1:10
    mu00= (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(lsh(1,xh0_2,xw0_2,10,8,lambda),0,ich0,icw0)-hdC(mu00,1)))^(crra); %l*(1/(w1_d(0.26,0)))*(1/(0.6856)^crrat) ; %0.45; %0.3885; 
    mu000=  (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(0,lsw(1,xh0_3,xw0_3,10,8,lambda),ich0,icw0)-hdC(mu000,1)))^(crra); %l*(1/(w1_d(0.26,0)))*(1/(0.6856)^crrat) ; %0.45; %0.3885; 
    mu0=  (lambda*ceh)*((1+((1-lambda)/lambda)^(1/crra))/(Yc(lsh(1,xh0_1,xw0_1,10,8,lambda),lsw(1,xh0_1,xw0_1,10,8,lambda),ich0,icw0)-hdC(mu0,1)))^(crra); %l*(1/(w1_d(0.258,0)))*(1/(0.69)^crrat); %0.42; %0.3902;0.6856
end

if (EQS.('inputs')==0)  
    fprintf('Creating own backup inputs')
        
    

    inputs0=zeros(T,T,I,I,I,3,3);
    for i=1:I
        for jh=1:I
            for jw=1:I 
                for th=1:T
                    for tw=1:T

                         inputs0(th,tw,jh,jw,i,1,:)=[mu00,xh0_2*100,xw0_2*100];
                         inputs0(th,tw,jh,jw,i,2,:)=[mu000,xh0_3*100,xw0_3*100];
                         inputs0(th,tw,jh,jw,i,3,:)=[mu0,xh0_1*100,xw0_1*100];
                    end
                end
            end
        end
    end
    inputs0_=inputs0;
        
else
    inputs0_=EQS.('inputs'); % substitutes 
end

     
if isempty(IN)        
    inputs0=inputs0_;
else
    inputs0=IN.('inputs');
    
end
inputs=inputs0;

for i=1:I
    for jh=1:I
        for jw=1:I 
            for th=1:T
                for tw=1:T
                    
if (hs_wfh==1) && ((th>T0)&& (tw<=T0))  && ((th<=T0)&& (tw>T0))   % these do not exist under household assignment - save time
    skip=1;
else    
                p=pp(i);
                
                dh=D(i,jh);
                dw=D(i,jw);
                if th>T0
                    dh=wfh_share*dh; %0;
                    [~,~,~,~,~,lowh]=matchdist(i,jh,th-T0,0,0,0,0,typeic,D,mm,JLs,betah);
                else
                    [~,~,~,~,~,lowh]=matchdist(i,jh,th,0,0,0,0,typeic,D,mm,JLs,betah);                    
                end
                if tw>T0
                    dw=wfh_share*dw; %0;
                    [~,~,~,~,~,loww]=matchdist(i,jw,tw-T0,0,0,0,0,typeic,D,mm,JLs, betah);
                else
                    [~,~,~,~,~,loww]=matchdist(i,jw,tw,0,0,0,0,typeic,D,mm,JLs, betah);
                end
                ich=1-lowh;
                icw=1-loww;
    
                rng(357);
                if piw_(dh,0)>0 && piw_(dh,0)<1 
                    fn=@(mu,xh,xw) [lsah_eq(mu,xh,xw,p,dh,0,ich,lambda) ,multC_eq(Yc(lsh(mu,xh,xw,dh,0,lambda),0,ich,icw),p,mu,xh,xw,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output1,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,1,:),options) ;
                end
                if piw_(dh,0)==1
                    %fprintf('xh=0')
                    fn=@(mu,Lh,xw) [lsah_eq_xh0(mu,Lh,xw,p,dh,0,ich,lambda) ,multC_eq(Yc(1-betah*dh - Lh,0,ich,icw),p,mu,0,xw,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output1,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,1,:),options) ;
                end
                if piw_(dh,0)==0
                    %fprintf('xw=0')
                    fn=@(mu,xh,xw) [lsah_eq_xw0(mu,xh,xw,p,dh,0,ich,lambda) ,multC_eq(Yc(lsh(mu,xh,xw,dh,0,lambda),0,ich,icw),p,mu,xh,xw,lambda)]; % one equation is just xw=0
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output1,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,1,:),options) ;
                end
                if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                   in_=inputs0(th,tw,jh,jw,i,1,:);
                    
                    if (i>1) || (jh>1) || (jw>1) || (th>1) || (tw>1)
                        [output1,~,EXITFLAG]=fsolve(fn, output1__,options);
                    else 
                        if VERBOSE
                            fprintf('Warning w0');    
                        end
                        [output1,~,EXITFLAG]=fsolve(fn,[in_(1),in_(2)*1.5,in_(3)],options0);  
                    end

                    if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                        if VERBOSE
                            fprintf('Warning w0 2');    
                        end
                        [output1,~,EXITFLAG]=fsolve(fn,[in_(1),in_(2)*2,in_(3)],options0);  
                        
                        if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                            
                            if VERBOSE
                                fprintf('Warning w0 3');    
                            end
                            [output1,~,EXITFLAG]=fsolve(fn,[mu00,in_(2)*1.1,in_(3)],options0);
                            if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))                            
                                if VERBOSE
                                    fprintf('Warning w0 3');    
                                end 
                                [output1,~,EXITFLAG]=fsolve(fn,[mu00,in_(2)*0.9,in_(3)],options0);
                                if ~isreal(output1) | output1(1)<=0 | output1(2)<=0 | lsh(output1(1),output1(2)/100,output1(3)/100,dh,0,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                    fprintf('labors: Warning w0')
                                    WARNINGS=WARNINGS+1;
                                    OUTC=999;
                                    inputs(th,tw,jh,jw,i,1,:)=inputs0_(th,tw,jh,jw,i,1,:);
                                    IN.('inputs')=real(inputs);
                                    return
                                end
                            end
                        end
                    end
                end
                 
                output1__=output1;
                inputs(th,tw,jh,jw,i,1,:)=output1;
                
                rng(357);
                if piw_(0,dw)>0 && piw_(0,dw)<1 
                    fn=@(mu,xh,xw) [lsaw_eq(mu,xh,xw,p,0,dw,icw,lambda) ,multC_eq(Yc(0,lsw(mu,xh,xw,0,dw,lambda),ich,icw),p,mu,xh,xw,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output2,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,2,:),options); 
                end
                if piw_(0,dw)==1
                    %fprintf('xh=0')
                    fn=@(mu,xh,xw) [lsaw_eq_xh0(mu,xh,xw,p,0,dw,icw,lambda) ,multC_eq(Yc(0,lsw(mu,xh,xw,0,dw,lambda),ich,icw),p,mu,xh,xw,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output2,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,2,:),options); 
                end
                if piw_(0,dw)==0
                    %fprintf('xw=0')
                    fn=@(mu,xh,Lw) [lsaw_eq_xw0(mu,xh,Lw,p,0,dw,icw,lambda) ,multC_eq(Yc(0,1-betah*dw - Lw,ich,icw),p,mu,xh,0,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output2,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,2,:),options); 
                end
                if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                    in_=inputs0(th,tw,jh,jw,i,2,:);
                    if (i>1) || (jh>1) || (jw>1) || (th>1) || (tw>1)
                        [output2,~,EXITFLAG]=fsolve(fn, output2__,options0);
                    else 
                        if VERBOSE
                            fprintf('Warning 0w');    
                        end
                        [output2,~,EXITFLAG]=fsolve(fn,[in_(1),in_(2),in_(3)*2],options0);
                    end

                    if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                        if VERBOSE
                            fprintf('Warning 0w 2');    
                        end
                        [output2,~,EXITFLAG]=fsolve(fn,[in_(1),in_(2),in_(3)*1.5],options0);
                    
                        if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))                       
                            if VERBOSE
                                fprintf('Warning 0w 3');    
                            end 
                            [output2,~,EXITFLAG]=fsolve(fn,[in_(1),in_(2),in_(3)*1.1],options0);
                            if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))                       
                                if VERBOSE
                                    fprintf('Warning 0w 4');    
                                end 
                                [output2,~,EXITFLAG]=fsolve(fn,[mu000,in_(2),in_(3)*0.9],options1);

                                if ~isreal(output2) | output2(1)<=0 | output2(2)<=0 | lsw(output2(1),output2(2)/100,output2(3)/100,0,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                    fprintf('labors: Warning 0w')
                                    WARNINGS=WARNINGS+1;
                                    OUTC=999;
                                    inputs(th,tw,jh,jw,i,2,:)=inputs0_(th,tw,jh,jw,i,2,:);
                                    IN.('inputs')=real(inputs);
                                    return
                                end
                            end
                        end
                    end
                end
                
                
                output2__=output2;
                inputs(th,tw,jh,jw,i,2,:)=output2;  
                
                rng(357);
                if piw_(dh,dw)>0 && piw_(dh,dw)<1
                    fn=@(mu,xh,xw) [lsb_eq(mu,xh,xw,p,dh,dw,ich,icw,lambda),multC_eq(Yc(lsh(mu,xh,xw,dh,dw,lambda),lsw(mu,xh,xw,dh,dw,lambda),ich,icw),p,mu,xh,xw,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output3,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,3,:),options);
                end
                if piw_(dh,dw)==1
                    %fprintf('xh=0')
                    fn=@(mu,Lh,xw) [lsb_eq_xh0(mu,Lh,xw,p,dh,dw,ich,icw,lambda),multC_eq(Yc(1-Lh-betah*dh,lsw(mu,0,xw,dh,dw,lambda),ich,icw),p,mu,0,xw,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output3,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,3,:),options);
                end
                if piw_(dh,dw)==0
                    %fprintf('xw=0')
                    fn=@(mu,xh,Lw) [lsb_eq_xw0(mu,xh,Lw,p,dh,dw,ich,icw,lambda),multC_eq(Yc(lsh(mu,xh,0,dh,dw,lambda),1-betah*dw - Lw,ich,icw),p,mu,xh,0,lambda)];
                    fn=@(x) fn(x(1),x(2)/100,x(3)/100);
                    [output3,~,EXITFLAG]=fsolve(fn,inputs0(th,tw,jh,jw,i,3,:),options);
                end

if ((jh==3) && (jw==2) && (th==1) && (tw==2) && (i==3))

    k=1;

end

                if (~isreal(output3)) | (output3(1)<=0) | (output3(2)<=0) |  (lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0)  | (lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0)  | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                    
                    
                    if (i>1) || (jh>1) || (jw>1) || (th>1) || (tw>1)
                        [output3,~,EXITFLAG]=fsolve(fn, output3__,options);
                    end 

                    if (~isreal(output3)) | (output3(1)<=0) | (output3(2)<=0) |  (lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0)  | (lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0)  | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                        if VERBOSE
                            fprintf('Warning');
                        end
                        in_= inputs0(th,tw,jh,jw,i,3,:);
                        if lsw(in_(1),in_(2)/100,in_(3)/100,dh,dw,lambda)<=0 
                            [output3,~,EXITFLAG]=fsolve(fn,[mu0,inputs0(th,tw,jh,jw,i,3,2)*0.1,inputs0(th,tw,jh,jw,i,3,3)],options0);
                        elseif lsh(in_(1),in_(2)/100,in_(3)/100,dh,dw,lambda)<=0
                            [output3,~,EXITFLAG]=fsolve(fn,[mu0,inputs0(th,tw,jh,jw,i,3,2),inputs0(th,tw,jh,jw,i,3,3)*0.1],options0);
                        else
                            [output3,~,EXITFLAG]=fsolve(fn,[mu0,reshape(inputs0_(th,tw,jh,jw,i,3,2:3),1,2)*1.1],options0);
                        end

                        if ~isreal(output3) | output3(1)<=0 | output3(2)<=0 |  lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0  | lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                            if VERBOSE
                                fprintf('Warning 2');
                            end

                            if lsw(in_(1),in_(2)/100,in_(3)/100,dh,dw,lambda)<=0 
                                [output3,~,EXITFLAG]=fsolve(fn,[mu0,inputs0(th,tw,jh,jw,i,3,2),inputs0(th,tw,jh,jw,i,3,3)*0.1],options1);
                            elseif lsh(in_(1),in_(2)/100,in_(3)/100,dh,dw,lambda)<=0
                                [output3,~,EXITFLAG]=fsolve(fn,[mu0,inputs0(th,tw,jh,jw,i,3,2)*0.1,inputs0(th,tw,jh,jw,i,3,3)],options1);
                            else
                                [output3,~,EXITFLAG]=fsolve(fn,[mu0,reshape(inputs0_(th,tw,jh,jw,i,3,2:3),1,2)*0.9],options1);
                            end

                            if ~isreal(output3) | output3(1)<=0 | output3(2)<=0 |  lsh(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0  | lsw(output3(1),output3(2)/100,output3(3)/100,dh,dw,lambda) <=0 | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                fprintf('labors: Warning2 ww')
                                WARNINGS=WARNINGS+1;
                                OUTC=999;
                                inputs(th,tw,jh,jw,i,3,:)=inputs0_(th,tw,jh,jw,i,3,:);
                                IN.('inputs')=real(inputs);
                                return
                            end
                        end

                    end
                end 
                
              output3__=output3;
              inputs(th,tw,jh,jw,i,3,:)=output3; % save the precise outputs aslater inputs
              
              %ROUNDING for outputs
              output1=round(output1,6); % to eliminate silly variation in outcomes with different initial guesses but same params
              output2=round(output2,8);
              output3=round(output3,8);
              
              
              output1(2:3)=output1(2:3)/100;
              output2(2:3)=output2(2:3)/100;
              output3(2:3)=output3(2:3)/100;
              Hc(th,tw,jh,jw,i,1)=hdC(output1(1),p);
              Hc(th,tw,jh,jw,i,2)=hdC(output2(1),p); 
              Hc(th,tw,jh,jw,i,3)=hdC(output3(1),p);
               
               vc(th,tw,jh,jw,i,1)=V(1,0,dh,0,AC(i),ch_fun(output1(1),lambda),cw_fun(output1(1),lambda),hdC(output1(1),p),lsh(output1(1),output1(2),output1(3),dh,0,lambda),0,...
                   xh_fun(output1(1),output1(2),output1(3),dh,0),xw_fun(output1(1),output1(2),output1(3),dh,0),0,ich,0,lambda); 
               vc(th,tw,jh,jw,i,2)=V(0,1,0,dw,AC(i),ch_fun(output2(1),lambda),cw_fun(output2(1),lambda),hdC(output2(1),p),0,lsw(output2(1),output2(2),output2(3),0,dw,lambda),...
                   xh_fun(output2(1),output2(2),output2(3),0,dw),xw_fun(output2(1),output2(2),output2(3),0,dw),0,0,icw,lambda); 
               vc(th,tw,jh,jw,i,3)=V(1,1,dh,dw,AC(i),ch_fun(output3(1),lambda),cw_fun(output3(1),lambda),hdC(output3(1),p),lsh(output3(1),output3(2),output3(3),dh,dw,lambda),lsw(output3(1),output3(2),output3(3),dh,dw,lambda),...
                   xh_fun(output3(1),output3(2),output3(3),dh,dw),xw_fun(output3(1),output3(2),output3(3),dh,dw),0,ich,icw,lambda); 
               
               
               if alloutput
                    lh(th,tw,jh,jw,i,1)=lsh(output1(1),output1(2),output1(3),dh,0,lambda) ; %output(1)*(output(1)>0)*(output(2)>0)*check; 
                    lw(th,tw,jh,jw,i,2)=lsw(output2(1),output2(2),output2(3),0,dw,lambda) ; %output(2)*(output(1)>0)*(output(2)>0)*check; 
                    lh(th,tw,jh,jw,i,3)=lsh(output3(1),output3(2),output3(3),dh,dw,lambda) ; %output(1)*(output(1)>0)*(output(2)>0)*check; 
                    lw(th,tw,jh,jw,i,3)=lsw(output3(1),output3(2),output3(3),dh,dw,lambda) ; %output(2)*(output(1)>0)*(output(2)>0)*check;
                    
                    mus(th,tw,jh,jw,i,1)=output1(1);
                    mus(th,tw,jh,jw,i,2)=output2(1);
                    mus(th,tw,jh,jw,i,3)=output3(1);
                    xh(th,tw,jh,jw,i,1)=xh_fun(output1(1),output1(2),output1(3),dh,0); 
                    xh(th,tw,jh,jw,i,2)=xh_fun(output2(1),output2(2),output2(3),0,dw); 
                    xh(th,tw,jh,jw,i,3)=xh_fun(output3(1),output3(2),output3(3),dh,dw); 
                    xw(th,tw,jh,jw,i,1)=xw_fun(output1(1),output1(2),output1(3),dh,0); 
                    xw(th,tw,jh,jw,i,2)=xw_fun(output2(1),output2(2),output2(3),0,dw); 
                    xw(th,tw,jh,jw,i,3)=xw_fun(output3(1),output3(2),output3(3),dh,dw); 
                    
                    ch(th,tw,jh,jw,i,1)=ch_fun(output1(1),lambda);
                    ch(th,tw,jh,jw,i,2)=ch_fun(output2(1),lambda);
                    ch(th,tw,jh,jw,i,3)=ch_fun(output3(1),lambda);
                    cw(th,tw,jh,jw,i,1)=cw_fun(output1(1),lambda);
                    cw(th,tw,jh,jw,i,2)=cw_fun(output2(1),lambda);
                    cw(th,tw,jh,jw,i,3)=cw_fun(output3(1),lambda);
                    uch(th,tw,jh,jw,i,1)=utilc(ch_fun(output1(1),lambda));
                    uch(th,tw,jh,jw,i,2)=utilc(ch_fun(output2(1),lambda));
                    uch(th,tw,jh,jw,i,3)=utilc(ch_fun(output3(1),lambda));
                    ucw(th,tw,jh,jw,i,1)=utilc(cw_fun(output1(1),lambda));
                    ucw(th,tw,jh,jw,i,2)=utilc(cw_fun(output2(1),lambda));
                    ucw(th,tw,jh,jw,i,3)=utilc(cw_fun(output3(1),lambda));
                    
                    ulh(th,tw,jh,jw,i,1)=utill(lh(th,tw,jh,jw,i,1),output1(2),dh);
                    ulh(th,tw,jh,jw,i,2)=utill(0,output2(2),0);
                    ulh(th,tw,jh,jw,i,3)=utill(lh(th,tw,jh,jw,i,3),output3(2),dh);
                    ulw(th,tw,jh,jw,i,1)=utill(0,output1(3),0);
                    ulw(th,tw,jh,jw,i,2)=utill(lw(th,tw,jh,jw,i,2),output2(3),dw);
                    ulw(th,tw,jh,jw,i,3)=utill(lw(th,tw,jh,jw,i,3),output3(3),dw);
                    
                    uh(th,tw,jh,jw,i,1)=utilh(1,dh,0,AC(i),ch_fun(output1(1),lambda),hdC(output1(1),p),lh(th,tw,jh,jw,i,1),...
                       xh_fun(output1(1),output1(2),output1(3),dh,0),xw_fun(output1(1),output1(2),output1(3),dh,0),0,ich); 
                    uh(th,tw,jh,jw,i,2)=utilh(0,0,dw,AC(i),ch_fun(output2(1),lambda),hdC(output2(1),p),0,...
                        xh_fun(output2(1),output2(2),output2(3),0,dw),xw_fun(output2(1),output2(2),output2(3),0,dw),0,ich); 
                    uh(th,tw,jh,jw,i,3)=utilh(1,dh,dw,AC(i),ch_fun(output3(1),lambda),hdC(output3(1),p),lh(th,tw,jh,jw,i,3),...
                        xh_fun(output3(1),output3(2),output3(3),dh,dw),xw_fun(output3(1),output3(2),output3(3),dh,dw),0,ich); 
                    uw(th,tw,jh,jw,i,1)=utilw(0,dh,0,AC(i),cw_fun(output1(1),lambda),hdC(output1(1),p),0,...
                        xh_fun(output1(1),output1(2),output1(3),dh,0),xw_fun(output1(1),output1(2),output1(3),dh,0),0,icw); 
                    uw(th,tw,jh,jw,i,2)=utilw(1,0,dw,AC(i),cw_fun(output2(1),lambda),hdC(output2(1),p), lw(th,tw,jh,jw,i,2),...
                        xh_fun(output2(1),output2(2),output2(3),0,dw),xw_fun(output2(1),output2(2),output2(3),0,dw),0,icw); 
                    uw(th,tw,jh,jw,i,3)=utilw(1,dh,dw,AC(i),cw_fun(output3(1),lambda),hdC(output3(1),p),lw(th,tw,jh,jw,i,3),...
                        xh_fun(output3(1),output3(2),output3(3),dh,dw),xw_fun(output3(1),output3(2),output3(3),dh,dw),0,icw); 

                    Y(th,tw,jh,jw,i,1)=Yc(lh(th,tw,jh,jw,i,1),0,ich,0);
                    Y(th,tw,jh,jw,i,2)=Yc(0,lw(th,tw,jh,jw,i,2),0,icw);
                    Y(th,tw,jh,jw,i,3)=Yc(lh(th,tw,jh,jw,i,3),lw(th,tw,jh,jw,i,3),ich,icw);
                    ux(th,tw,jh,jw,i,1)=utilx(xh_fun(output1(1),output1(2),output1(3),dh,0),xw_fun(output1(1),output1(2),output1(3),dh,0),0,dh,0);
                    ux(th,tw,jh,jw,i,2)=utilx(xh_fun(output2(1),output2(2),output2(3),0,dw),xw_fun(output2(1),output2(2),output2(3),0,dw),0,0,dw);
                    ux(th,tw,jh,jw,i,3)=utilx(xh_fun(output3(1),output3(2),output3(3),dh,dw),xw_fun(output3(1),output3(2),output3(3),dh,dw),0,dh,dw);
                    ucl(th,tw,jh,jw,i,1)=utilcl(dh,0);
                    ucl(th,tw,jh,jw,i,2)=utilcl(0,dw);
                    ucl(th,tw,jh,jw,i,3)=utilcl(dh,dw);
                    uho(th,tw,jh,jw,i,1)=utilho(Hc(th,tw,jh,jw,i,1));
                    uho(th,tw,jh,jw,i,2)=utilho(Hc(th,tw,jh,jw,i,2));
                    uho(th,tw,jh,jw,i,3)=utilho(Hc(th,tw,jh,jw,i,3));
                    
                    uXih(th,tw,jh,jw,i,1)=uXi(lh(th,tw,jh,jw,i,1),ich);
                    uXih(th,tw,jh,jw,i,2)=uXi(0,0);
                    uXih(th,tw,jh,jw,i,3)=uXi(lh(th,tw,jh,jw,i,3),ich);
                    uXiw(th,tw,jh,jw,i,1)=uXi(0,0);
                    uXiw(th,tw,jh,jw,i,2)=uXi(lw(th,tw,jh,jw,i,2),icw);
                    uXiw(th,tw,jh,jw,i,3)=uXi(lw(th,tw,jh,jw,i,3),icw);
                    ich_(th,tw,jh,jw,i,1)=ich;
                    ich_(th,tw,jh,jw,i,2)=ich;
                    ich_(th,tw,jh,jw,i,3)=ich;
                    icw_(th,tw,jh,jw,i,1)=icw;
                    icw_(th,tw,jh,jw,i,2)=icw;
                    icw_(th,tw,jh,jw,i,3)=icw;
                    
                    piw_dd(th,tw,jh,jw,i,1)=piw_(dh,0);
                    piw_dd(th,tw,jh,jw,i,2)=piw_(0,dw);
                    piw_dd(th,tw,jh,jw,i,3)=piw_(dh,dw);
               end
                
end               

                end
            end
        end
    end
end


OUTC.('vc') = vc-vc(1,1,1,1,1,1); % i.e. do not use for any welfare comparisons, just for generating loc choices!!!

OUTC.('Hc') = Hc;
OUTC.('lambda') = lambda;
OUTC.('p') = pp;
IN.('inputs')=real(inputs);

if alloutput
    OUTC.('vc_scale') = vc;
    OUTC.('inputs0')=inputs0; 
    OUTC.('inputs')=real(inputs);
    OUTC.('Yc') = Y;
    OUTC.('xh') = xh;
    OUTC.('lh') = lh;
    %OUTC.('ich') = ich; % nonsense
    OUTC.('xw') = xw;
    OUTC.('lw') = lw;
    %OUTC.('icw') = icw;
    OUTC.('mus') = mus;
    OUTC.('ch') = ch;
    OUTC.('cw') = cw;
    OUTC.('uh') = uh;
    OUTC.('uw') = uw;
    OUTC.('uch') = uch;
    OUTC.('ucw') = ucw;
    OUTC.('ulh') = ulh;
    OUTC.('ulw') = ulw;
    OUTC.('ux') = ux;
    OUTC.('ucl') = ucl;
    OUTC.('uho') = uho;
    OUTC.('uXih') = uXih;
    OUTC.('uXiw') = uXiw;
    OUTC.('ich_') = ich_;
    OUTC.('icw_') = icw_;
    OUTC.('epsC') = 0;
end

end