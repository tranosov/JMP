function OUTS=laborsS(pp,alloutput,EQS,PARREST)
%global EQS JLs D AS typeic WARNINGS mm params
global WARNINGS
global INS
rng(357)

typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
AS=PARREST.('AS');
betah=PARREST.('betah');
crra=PARREST.('crra');
ces=PARREST.('ces');

if pp==99
    params=PARREST.('params');
    pp=[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ];  % what was this for?
end

lssh_eq=EQS.('lssh_eq');  
multS_eq=EQS.('multS_eq');
lssh=EQS.('lssh') ;
Ys=EQS.('Ys');
us=EQS.('us');
cs=EQS.('cs');
hdS=EQS.('hdS');

N=size(JLs);
I=N(1);
T=N(2);
vs=zeros(T,I,I);
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
mu0=ces*(1/( Ys(0.23,0)-1))^(crra); 
    
if isempty(INS)  
    if EQS.('inputsS')==0        
        inputs0=zeros(T,I,I,2);
        for j=1:I
            for i=1:I
                for t=1:T
                    inputs0(t,j,i,:)=[mu0,x0];

                end
            end
        end

    else 
        inputs0=EQS.('inputsS');
    end
else
    inputs0=INS.('inputsS'); % within estimation use last inputs
end

TOL=10^(-15);
options = optimoptions('fsolve','MaxIter',5000,'MaxFunctionEvaluations',5000,...
        'FunctionTolerance',TOL,'Display','off','Algorithm','levenberg-marquardt','ScaleProblem','jacobian');
inputs=zeros(T,I,I,2);
for j=1:I
    for i=1:I
        for t=1:T
            d=D(i,j);
            p=pp(i);
            [~,~,~,~,~,low]=matchdist(i,j,t,0,0,0,0,typeic,D,mm,JLs,betah);
            ic=1-low;
            ics(t,j,i)=ic;
            fn=@(mu,x) [lssh_eq(mu,x,p,d,ic) ,multS_eq(Ys(lssh(mu,x,d),ic),p,mu,x)];
            fn=@(in) fn(in(1),in(2));
            [output,~,EXITFLAG]=fsolve(fn,inputs0(t,j,i,:),options); 
            if (~isreal(output) ) | (output(1)<=0 )| ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                    %fprintf('Warning')
                    %inputs0(t,j,i,:).*sign(fn(inputs0(t,j,i,:))).*[-1,1]*1.1
                    [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,inputs0(t,j,i,:).*sign(fn(inputs0(t,j,i,:))).*[-1,1]*1.1,options);
                    if (~isreal(output)) | (output(1)<=0) | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                            [output,FVAL,EXITFLAG,OUTPUT]=fsolve(fn,[mu0,x0],options);
                            if (~isreal(output)) | (output(1)<=0) | ((EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3)&& (EXITFLAG~=4))
                                fprintf('Warning on singles - 2')
                                WARNINGS=WARNINGS+1;
                                
                            end      
                    end
            end
            
            %fn=@(ls) lssw_eq(ls,p,d,ic);
            %lssw(t,j,i)=fsolve(fn,ls0,options); 
            %xsw(t,j,i)=xsw_fun(lssw(t,j,i),d);
            inputs(t,j,i,:)=output;
            Hs(t,j,i)=hdS(output(1),p);
            vs(t,j,i)= us(1,d,AS(i),cs(output(1)),Hs(t,j,i),lssh(output(1),output(2),d),output(2),0,ic);
            
            if alloutput
                mus(t,j,i)=output(1);
                xs(t,j,i)=output(2);
                ics(t,j,i)=ic;
                lss(t,j,i)=lssh(output(1),output(2),d);
                ucs(t,j,i)=utilcs(cs(output(1)));
                conss(t,j,i)=cs(output(1));
                uhs(t,j,i)=utilhs(Hs(t,j,i));
                uls(t,j,i)=utills(lssh(output(1),output(2),d),output(2),d);
                ucls(t,j,i)=utilcls(d);
                uxs(t,j,i)=utilxs(output(2),0);
                Yss(t,j,i)=Ys(lssh(output(1),output(2),d),ic);
            end
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
   
end

end