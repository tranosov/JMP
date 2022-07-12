
function [output,EXITFLAG,time]=solve(x0,EQS,PARREST)
global WARNINGS RESC VERBOSE
time=0;
rng(357)
if VERBOSE
    iter_="iter";
else
    iter_="off";
end

if isempty(RESC)
    RESC=10^6;
end

%{
WARNINGS=0;
    kk=1
    cl=F(x0) % in lss: this is already off
    while (sum(abs(cl(end)))>1) && (kk<5)
        x0=x0 +cl.*([0,0,0,(RESC)*PARREST.('sigmam')/(1000*kk)]) 
        cl=F(x0) 
        kk=kk+1;
    end
    kk=1;
    while (sum(abs(cl(1:end-1)))>1) && (kk<5)
        A=1/(1000*kk);
        x0=x0 +cl.*(A*[1,1,1,0])
        [cl]=F(x0) 
        kk=kk+1;
    end  
    WARNINGS=0;

toc
time=toc;

%}      
         
tic
global mup
if isempty(mup)
    x0_=x0(end);
else
    if mup==1
        x0_=[max(0.001,x0(end)*0.95),x0(end)*1.01]; % more men in mm, clm(x0(end)) should be<0, lambda down to incite women
        x02_=[max(0.001,x0(end)*0.5),x0(end)*1.01];
        x03_=[max(0.001,x0(end)*0.5),x0(end)*1.1];
    end
    if mup==(-1)
        x0_=[0.99*x0(end),min(x0(end)*1.05,0.999)]; % less men in mm, clm(x0(end)) should be>0, lambda up to incite men
        x02_=[0.99*x0(end),min(x0(end)*2,0.999)];
        x03_=[0.9*x0(end),min(x0(end)*2,0.999)];
    end
end 

Fm=@(x) Clearing_justmm(x,EQS,PARREST); 
tol=10^(-3); % stricter than overall
if norm(Fm(x0(end)))^2 >tol 
    optionsz = optimset('TolX',tol,'Display',iter_);
    try
        out=fzero(Fm,x0_,optionsz); % something is off here...
    catch
        try
            out=fzero(Fm,x02_,optionsz); 
        catch
             try
                fprintf('unexpected behavior at solving for lambda in solve')
                out=fzero(Fm,x03_,optionsz); 
            catch
                fprintf('Failed to resolve for lambda in solve');
                output=[999,999,999,999];
                EXITFLAG=999;
                time=time+toc;
                return
            end
        end
    end
    x0(end)=out;
end 

params=PARREST.('params');
LA0_=params('LA0','value');
params('LA0','value')={x0(end)/(RESC)};
PARREST.('params')=params;
time=time+toc;
if VERBOSE
toc
end

tic
Fp=@(x) Clearing(x,EQS,PARREST); 
tol=1;
options = optimoptions('fsolve','MaxIter',500,'MaxFunctionEvaluations',5000,...
           'FunctionTolerance',tol*10^(-1), 'Display',iter_,...
           'Algorithm','trust-region'); 
if norm(Fp(x0(1:end-1)))^2 >tol*10^(-1)   
    x0(1:end-1)=fsolve(Fp,x0(1:end-1),options);
end
params('LA0','value')=LA0_;
PARREST.('params')=params;
time=time+toc;
if VERBOSE
toc
end

            
WARNINGS=0;
F=@(x) Clearing_withmm(x,EQS,PARREST); 
cl=F(x0); 
if WARNINGS>0
    output=[999,999,999,999];
    EXITFLAG=999;
    fprintf('ISSUES AT EVALUATION THE CONTINUOUS VARIABLES in solve.');

else 
    if norm(cl)^2 >tol  
                tic
if VERBOSE
fprintf('Have to do joint solution...')
end
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0,options);
                % does not really work - check whether I can help it somehow?

                % first solve lambda and only then solve prices? the feedback
                % seems to be mostly from lambda to prices and not the other
                % way around
                if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)
                    fprintf('Lowering standards in solve...');
                    kkk=1;
                    while (sum(abs(cl))>0.5) && (kkk<3)
                        kk=1;
                        while (sum(abs(cl(end)))>0.1) && (kk<10)
                            x0=x0 +cl.*([0,0,0,(RESC)*PARREST.('sigmam')/(300*kk)]);
                            cl=F(x0);
                            if VERBOSE
                                cl
                            end
                            kk=kk+1;
                        end
                        kk=1;
                        while (sum(abs(cl(1:end-1)))>0.4) && (kk<10)
                            A=1/(1000*kk);
                            x0=x0 +cl.*(A*[1,1,1,0]);
                            [cl]=F(x0); 
                            if VERBOSE
                                cl
                            end
                            kk=kk+1;
                        end  
                        WARNINGS=0;
                        kkk=kkk+1;
                    end
                    options = optimoptions('fsolve','MaxIter',500,'Display',iter_,'MaxFunctionEvaluations',5000,...
                    'FunctionTolerance',10^(0), 'Algorithm', 'trust-region');
                    [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0,options);

                end
                if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)

                    WARNINGS=WARNINGS+0.5; % had to lower standards
                    fprintf('Trying again');
                    [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0*0.9,options);

                end
                if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)

                    fprintf('Trying again again');
                    [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0*1.1,options);

                end
                if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)

                    fprintf('Trying again again again');
                    [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,[1,1,1,0.5],options);

                    if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)              
                        WARNINGS=WARNINGS+1;
                    end

                end
                
                if VERBOSE
                toc
                end

                time=time+toc;
    else
        output=x0;
        EXITFLAG=1;
    end

end
        
        
end
        