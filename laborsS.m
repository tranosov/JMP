function OUTS=laborsS(pp,alloutput,EQS,PARREST)
%global EQS JLs D AS typeic WARNINGS mm params
global WARNINGS VERBOSE
global INS

typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
AS=PARREST.('AS');
betah=PARREST.('betah');
crra=PARREST.('crra');
ces=PARREST.('ces');
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

lssh_eq=EQS.('lssh_eq');  
multS_eq=EQS.('multS_eq');
lssh=EQS.('lssh') ;
Ys=EQS.('Ys');
us=EQS.('us');
cs=EQS.('cs');
hdS=EQS.('hdS');
Leissh=EQS.('Leissh');
mus0=EQS.('mus0');  
uXi=EQS.('uXi');

wfh=PARREST.('wfh');
T0=T;
if wfh>0
    T=T0*2;
end

vs=zeros(T,I,I,2);
Hs=zeros(T,I,I);
if alloutput
    utilhs=EQS.('utilhs') ;
    utilcs=EQS.('utilcs');
    utilxs=EQS.('utilxs') ;
    utills=EQS.('utills') ;
    utilcls=EQS.('utilcls') ;
    lss=zeros(T,I,I);
    ics=zeros(T,I,I);
    %lssw=zeros(T,I,I); - not implemented to be different
    xs=zeros(T,I,I);
    mus=zeros(T,I,I);
    ucs=zeros(T,I,I);
    conss=zeros(T,I,I);
    uhs=zeros(T,I,I);
    uls=zeros(T,I,I);
    ucls=zeros(T,I,I);
    uxs=zeros(T,I,I);
    Yss=zeros(T,I,I);
end

x0=0.05;
mu0=ces*(1/( Ys(lssh(1,x0,8),0)-1))^(crra);

if EQS.('inputsS')==0        
        inputs0_=zeros(T,I,I,2);
        for j=1:I
            for i=1:I
                for t=1:T
                    inputs0_(t,j,i,:)=[mu0,x0];

                end
            end
        end

    else 
        inputs0_=EQS.('inputsS');
end
if isempty(INS)  
    inputs0=inputs0_;
else
    inputs0=INS.('inputsS'); % within estimation use last inputs
end
inputs=inputs0;

TOL=10^(-15);
STEPTOL=10^(-10);
options = optimoptions('fsolve','MaxIter',500,'MaxFunctionEvaluations',500,...
        'FunctionTolerance',TOL,'Display','off','Algorithm','levenberg-marquardt','ScaleProblem','jacobian',...
        'StepTolerance', STEPTOL);

for j=1:I
    for i=1:I
        for t=1:T
            d=D(i,j);
            p=pp(i);

            if t>T0
                d=wfh_share*d; %0;
                [~,~,~,~,~,low]=matchdist(i,j,t-T0,0,0,0,0,typeic,D,mm,JLs,betah);
            else
                [~,~,~,~,~,low]=matchdist(i,j,t,0,0,0,0,typeic,D,mm,JLs,betah);
            end
            ic=1-low;
            ics(t,j,i)=ic;
            
            fn_=@(mu,x) [lssh_eq(mu,x,p,d,ic) ,multS_eq(Ys(lssh(mu,x,d),ic),p,mu,x)];
            fn=@(in) fn_(in(1),in(2));
            
            rng(357);
            [output,~,EXITFLAG]=fsolve(fn,inputs0(t,j,i,:),options); 
            output_=output;
            if (~isreal(output) ) || (output(1)<=0 )|| ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4)) || (numel(output)~=numel(inputs(t,j,i,:)))
                    if VERBOSE
                        fprintf('laborS: warning')
                    end

                    [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,inputs0(t,j,i,:).*sign(fn(inputs0(t,j,i,:))).*[-1,1]*2,options);
                    
                    if (~isreal(output)) || (output(1)<=0) || ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))|| (numel(output)~=numel(inputs(t,j,i,:)))
                               [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,inputs0(t,j,i,:).*sign(fn(inputs0(t,j,i,:))).*[-1,1]*0.5,options);
                               
                        if (~isreal(output)) || (output(1)<=0) || ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))|| (numel(output)~=numel(inputs(t,j,i,:)))
                            options0 = optimoptions('fsolve','MaxIter',5000,'MaxFunctionEvaluations',5000,...
                                'FunctionTolerance',TOL,'Display','off','Algorithm','trust-region','StepTolerance', STEPTOL);
                            [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,[mu0,x0],options0);
                            
                            if (~isreal(output)) | (output(1)<=0) | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4)) || (numel(output)~=numel(inputs(t,j,i,:)))
                                fprintf('laborsS:  warning 2')
                                OUTS=999;
                                inputs(t,j,i,:)=inputs0_(t,j,i,:);
                                INS.('inputsS')=real(inputs);
                                WARNINGS=WARNINGS+1;
                                return
                            end  
                        end
                        
                        
                    end
            end
            
            %fn=@(ls) lssw_eq(ls,p,d,ic);
            %lssw(t,j,i)=fsolve(fn,ls0,options); 
            %xsw(t,j,i)=xsw_fun(lssw(t,j,i),d);
            inputs(t,j,i,:)=output;
            Hs(t,j,i)=hdS(output(1),p);
            vs(t,j,i,2)= us(1,d,AS(i),cs(output(1)),Hs(t,j,i),lssh(output(1),output(2),d),output(2),0,ic);
            
            if alloutput
                mus(t,j,i)=output(1);
                phs(t,j,i)=Hs(t,j,i)*p;
                xs(t,j,i)=output(2);
                ics(t,j,i)=ic;
                lss(t,j,i)=lssh(output(1),output(2),d);
                ucs(t,j,i)=utilcs(cs(output(1)));
                Yss(t,j,i)=Ys(lssh(output(1),output(2),d),ic);
                conss(t,j,i)=cs(output(1));
                uhs(t,j,i)=utilhs(Hs(t,j,i));
                uls(t,j,i)=utills(lssh(output(1),output(2),d),output(2),d);
                ucls(t,j,i)=utilcls(d);
                uxs(t,j,i)=utilxs(output(2),0);

                uxi(t,j,i)=uXi(lssh(output(1),output(2),d),ic);
   
            end

        
        % what if not working?    
        %mu_=mus0(p);
        %fn=@(x) [Leissh(mu_,x)-(1-x)];
        %[x_,~,EXITFLAG]=fsolve(fn,x0,options); 

        %vs(t,j,i,1)= us(0,d,AS(i),cs(mu_),Hs(t,j,i),0,x_,0,ic);


        end
    end
end
OUTS = struct();
OUTS.('vs') = vs;
OUTS.('Hs') = Hs;
OUTS.('p') = pp;
INS.('inputsS')=real(inputs);

if alloutput
    OUTS.('inputs0')=inputs0; 
    OUTS.('inputs')=real(inputs);
    OUTS.('Ys') = Yss;
    OUTS.('xs') = xs;
    OUTS.('lss') = lss;
    OUTS.('ics') = ics;
    OUTS.('ucs') = ucs;
    OUTS.('conss') = conss;
    OUTS.('uls') = uls;
    OUTS.('uxs') = uxs;
    OUTS.('ucls') = ucls;
    OUTS.('uhs') = uhs; 
    OUTS.('mus') = mus; 
    OUTS.('phs') = phs; 
    OUTS.('epsS') = 0;
    OUTS.('uxi') = uxi;
end



end