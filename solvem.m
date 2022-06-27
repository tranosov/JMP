
function [output,EXITFLAG,time]=solvem(x0,EQS,PARREST)
global WARNINGS RESC VERBOSE
rng(357)

F=@(x) Clearing_justmm(x,EQS,PARREST); 

WARNINGS=0;
    cl=F(x0)
    kk=1;
    while (sum(abs(cl))>1) & (kk<5) %1
        x0=x0 +(RESC)*PARREST.('sigmam')/(1000*kk);
        cl=F(x0);
        if VERBOSE
            cl=cl
        end
        kk=kk+1;
        WARNINGS=0;
    end

toc
time=toc;

        
        
if WARNINGS>0
    fprintf('ISSUES AT EVALUATION THE CONTINUOUS VARIABLES in solvem.');
    output=[999];
    EXITFLAG=999;

else
            tic
            optionsz = optimset('TolX',tol*10^(-6),'Display','off');
            WARNINGS=0;
            if VERBOSE
                optionsz = optimset('TolX',tol*10^(-6),'Display','iter');
            end             
            [output,FVAL,EXITFLAG,OUTPUT]= fzero(F,x0,optionsz);
            
            
            if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)
                OUTPUT
                fprintf('Trying again');
                kk=1;
                    while (sum(abs(cl))>0.1) && (kk<10)
                        x0=x0 +cl.*(RESC)*PARREST.('sigmam')/(200*kk);
                        kk=kk+1;
                        cl=F(x0)
                    end
                if VERBOSE
                    options = optimoptions('fsolve','MaxIter',500,'Display','iter','MaxFunctionEvaluations',5000,...
                'FunctionTolerance',10^(-6), 'Display','iter','Algorithm', 'trust-region');
                else
                    options = optimoptions('fsolve','MaxIter',500,'Display','iter','MaxFunctionEvaluations',5000,...
                                'FunctionTolerance',10^(-6), 'Display','off','Algorithm', 'trust-region');
                end 
                
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0,options);
                EXITFLAG
            end

            if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)
                OUTPUT
                WARNINGS=WARNINGS+0.5; % had to lower standards
                fprintf('Trying again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0*0.9,options);
                EXITFLAG
            end
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)
                OUTPUT
                fprintf('Trying again again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0*1.1,options);
                EXITFLAG
            end
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)
                OUTPUT
                fprintf('Trying again again again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,[0.5],options);
                EXITFLAG
                if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)              
                    WARNINGS=WARNINGS+1;
                    fprintf('Solvem saved many times.')
                    EXITFLAG=999;
                end
            end
            toc;
            
            time=time+toc;
end
        
        
end
        