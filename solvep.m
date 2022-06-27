
function [output,EXITFLAG,time]=solvep(lp0,EQS,PARREST)
global WARNINGS VERBOSE
rng(357);
if VERBOSE
    iter_="iter";
else
    iter_="off";
end


F=@(x) Clearing(x,EQS,PARREST); 
params=PARREST.('params');
tic
WARNINGS=0;
    cl0=F(lp0);
    cl=cl0;
    kk=1;
    while (sum(abs(cl))>5) && (kk<5)
            A=500*table2array(params('crrah_','value'));
            if VERBOSE
                cl=cl
            end
            if WARNINGS>0
                lp0=lp0+cl0./(A*[1,1,1]);
            else
                lp0=lp0+cl./(A*[1,1,1]); % roughly correct initial guess
            end
            
        kk=kk+1;
        WARNINGS=0;
        cl=F(lp0);
    end

time=toc;
if VERBOSE
    toc
    cl=cl
end
  
        if WARNINGS>0
            fprintf('ISSUES AT EVALUATION THE CONTINUOUS VARIABLES in solvep.');
            output=[999,999,999];
            EXITFLAG=999;
            return
        else
            tic
            
            WARNINGS=0;
            options = optimoptions('fsolve','MaxIter',50,'MaxFunctionEvaluations',500,...
                       'FunctionTolerance',10^(-3), 'Display',iter_,...
                       'FunValCheck','on','Algorithm','trust-region'); %'trust-region-dogleg' 'Display','final'


            [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0,options);
            if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)

                %fprintf('Trying again');
                cl=F(output);
                kk=1;
                    while (sum(abs(cl))>0.5) && (kk<10)
                        A=1000*table2array(params('crrah_','value'));
                        lp0=lp0+cl./(A*[1,1,1]); % roughly correct initial guess
                        kk=kk+1;
                        cl=F(lp0);
                        if VERBOSE
                            kk
                            cl=cl
                        end
                    end

                options = optimoptions('fsolve','MaxIter',500,'Display',iter_,'MaxFunctionEvaluations',5000,...
                'FunctionTolerance',10^(0), 'Algorithm', 'trust-region');
                
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0,options);

            end
            if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)

                WARNINGS=WARNINGS+0.5; % had to lower standards
                %fprintf('Trying again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0*0.9,options);

            end
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)

                %fprintf('Trying again again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0*1.1,options);

            end
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)

                %fprintf('Trying again again again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,[1,1,1],options);

                if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)              
                    WARNINGS=WARNINGS+1;
                    fprintf('Solvep failed many times.')
                    EXITFLAG=999;
                end

            end

            time=time+toc;
            if VERBOSE
            toc
            end
            
            
        end
        
        
end
        