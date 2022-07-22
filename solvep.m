
function [output,EXITFLAG,time]=solvep(lp0,EQS,PARREST,tolmult)
global WARNINGS VERBOSE
rng(357);
if VERBOSE
    iter_="iter";
else
    iter_="off";
end

TOL=10^(-1)*tolmult; 
F=@(x) Clearing(x,EQS,PARREST); 
params=PARREST.('params');
tic
WARNINGS=0;
    cl0=F(lp0);
    cl=cl0;
    kk=1;
    I=size(lp0,2);
        if WARNINGS==0
            while (kk<15) && (sum(abs(cl))>TOL) 
                    A=800/table2array(params('crrah_','value'));
                    if VERBOSE
                        cl=cl
                    end
                    if WARNINGS>0
                        lp0=lp0+cl0./(A*ones(1,I)); % go back
                    else
                        lp0=lp0+cl./(A*ones(1,I)); % roughly correct initial guess
                    end

                kk=kk+1;
                WARNINGS=0;
                cl=F(lp0);
            end
        end

time=toc;
if VERBOSE
    toc
    cl=cl
end

%STEPTOL=10^(-8);

        if WARNINGS>0
            fprintf('ISSUES AT EVALUATION THE CONTINUOUS VARIABLES in solvep.');
            output=reshape(999, size(lp0));
            EXITFLAG=999;
            return
        elseif norm(cl)^2 <=TOL
            output=lp0;
            EXITFLAG=1;
        else
            tic
            
            WARNINGS=0;
            options = optimoptions('fsolve','MaxIter',50,'MaxFunctionEvaluations',500,...
                       'FunctionTolerance',TOL,   'Display',iter_,...
                       'FunValCheck','on','Algorithm','levenberg-marquardt'); %'trust-region-dogleg' 'Display','final'
%'StepTolerance', STEPTOL?
            try
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0,options);
            catch
                EXITFLAG=999;
                output=lp0;
            end

            if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)

                %fprintf('Trying again');
                cl=F(output);
                kk=1;
                    while (sum(abs(cl))>0.5) && (kk<10)
                        A=1000/table2array(params('crrah_','value'));
                        lp0=lp0+cl./(A*[1,1,1]); % roughly correct initial guess
                        kk=kk+1;
                        cl=F(lp0);
                        if VERBOSE
                            kk
                            cl=cl
                        end
                    end

                options = optimoptions('fsolve','MaxIter',500,'Display',iter_,'MaxFunctionEvaluations',5000,...
                'FunctionTolerance',   TOL, 'Algorithm', 'trust-region');
                if norm(cl)^2 <=TOL
                    output=lp0;
                    EXITFLAG=1;
                else
                    [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0,options);
                end

            end
            if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)

                WARNINGS=WARNINGS+0.5; % had to lower standards
                fprintf('Lowering standards in solvep');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0*0.9,options);

            end
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)              
                WARNINGS=WARNINGS+1;
                fprintf('Solvep failed many times.')
                FVAL=FVAL
                OUTPUT
                EXITFLAG=999;
            end
            
            
            %{
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)

                %fprintf('Trying again again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,lp0*1.1,options);

            end
            if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)

                %fprintf('Trying again again again');
                [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,ones(size(lp0)),options);

            

            end
%}
            time=time+toc;
            if VERBOSE
                toc
            end
            
            
        end
        
output=round(output,6);        
end
        