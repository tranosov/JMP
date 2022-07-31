
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

if mup==1
    x0_=[max(0.001*RESC,x0(end)*0.9),x0(end)*1.01]; % more men in mm, clm(x0(end)) should be<0, lambda down to incite women
    x02_=[max(0.001*RESC,x0(end)*0.5),x0(end)*1.05];
    %x03_=[max(0.001*RESC,x0(end)*0.5),x0(end)*1.1];
end
if mup==(-1)
    x0_=[0.99*x0(end),min(x0(end)*1.1,0.999*RESC)]; % less men in mm, clm(x0(end)) should be>0, lambda up to incite men
    x02_=[0.95*x0(end),min(x0(end)*1.5,0.999*RESC)];
    %x03_=[0.9*x0(end),min(x0(end)*1.5,0.999*RESC)];
end

Fm=@(x) Clearing_justmm(x,EQS,PARREST); 
tol=10^(-8); % stricter than overall
if norm(Fm(x0(end)))^2 >tol 
    optionsz = optimset('TolX',tol,'Display',iter_);
    if isempty(mup)
        try
           out=fzero(Fm,x0(end),optionsz); 
        catch
            fprintf('Failed to resolve for lambda in solve');
            output=reshape(999,size(x0));
            EXITFLAG=991;
            time=time+toc;
            return
        end
        
    else
        fmx2=Fm(x0_(2));
        fmx1=Fm(x0_(1));
        if fmx2*fmx1<0
            out=fzero(Fm,x0_,optionsz);
            fmout=Fm(out);
            if fmout<-0.01
                if fmout*fmx2<-10^(-3)
                    out=fzero(Fm,[out,x0_(2)],optionsz);
                elseif  fmout*fmx1<-10^(-3)
                    out=fzero(Fm,[out,x0_(1)],optionsz);
                else
                    fprintf('Fm(out) %d <0 , fmx2 %d>0 or fmx1 %d>0\n',fmout,fmx2,fmx1)
                    fprintf('Yet somehow I cannot solve\n')
                end
            elseif fmout>0.01
                if fmout*fmx1<-10^(-3)
                    out=fzero(Fm,[x0_(1),out],optionsz);
                elseif fmout*fmx2<-10^(-3)
                    out=fzero(Fm,[out,x0_(2)],optionsz);
                else
                    fprintf('Fm(out) %d >0 , fmx1 %d>0 or fmx2 %d>0 \n',fmout,fmx1,fmx2)
                    fprintf('Yet somehow I cannot solve\n')
                end
            end
        else
            optionsz = optimset('TolX',tol,'Display','iter');
            fprintf('Wider net on lambda?')
            fmx2=Fm(x02_(2));
            fmx1=Fm(x02_(1));
            if fmx2*fmx1<-(10^(-3))
                out=fzero(Fm,x02_,optionsz); 
                fmout=Fm(out);
                if fmout<-0.01
                    if fmout*fmx2<-10^(-3)
                        out=fzero(Fm,[out,x0_(2)],optionsz);
                    elseif fmout*fmx2<-10^(-3)
                        out=fzero(Fm,[out,x0_(2)],optionsz);
                    else
                        fprintf('Fm(out) %d >0 , fmx1 %d>0 or fmx2 %d>0 \n',fmout,fmx1,fmx2)
                        fprintf('Yet somehow I cannot solve\n')
                    end
                elseif fmout>0.01
                    if fmout*fmx1<-10^(-3)
                        out=fzero(Fm,[x0_(1),out],optionsz);
                    elseif fmout*fmx2<-10^(-3)
                        out=fzero(Fm,[out,x0_(2)],optionsz);
                    else
                        fprintf('Fm(out) %d >0 , fmx1 %d>0 or fmx2 %d>0 \n',fmout,fmx1,fmx2)
                        fprintf('Yet somehow I cannot solve\n')
                    end
                end
            else
                fprintf('Unexpected behavior for lambda in solve');
                try
                   out=fzero(Fm,x0(end),optionsz); 
                catch
                    fprintf('Failed to resolve for lambda in solve');
                    output=reshape(999,size(x0));
                    EXITFLAG=991;
                    time=time+toc;
                    return
                end
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
cl=Fp(x0(1:end-1));    
if norm(cl)^2 >10^(-1)*tol
    tolmult=10^(-1);
    [x0(1:end-1),EXITFLAG,time_]=solvep(x0(1:end-1),EQS,PARREST,tolmult);
    time=time+time_;
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
    output=reshape(999,size(x0));
    EXITFLAG=999;
    fprintf('ISSUES AT EVALUATION THE CONTINUOUS VARIABLES in solve.');
    return
else 
    
    if norm(cl)^2 >tol  
                tic
        %if VERBOSE
        fprintf('Have to do joint solution...')
        cl
        %end

    TOL=tol; 
    %STEPTOL=10^(-8);
    options = optimoptions('fsolve','MaxIter',500,'MaxFunctionEvaluations',500,...
                           'FunctionTolerance',TOL,   'Display','iter',...
                           'FunValCheck','on','Algorithm','trust-region'); % switch to fminsearch? not yet.
    
    
    
                %[output,FVAL,EXITFLAG,OUTPUT]= fminsearch(F_ob,x0,options);
                % does not really work - check whether I can help it somehow?

                % first solve lambda and only then solve prices? the feedback
                % seems to be mostly from lambda to prices and not the other
                % way around
                    kkk=1;
                    while (sum(abs(cl))>0.5) && (kkk<3)
                        kk=1;
                        fprintf('Helping lambda...');
                        while (sum(abs(cl(end)))>0.1) && (kk<15)
                            x0=x0 +cl.*([zeros(size(x0(1:end-1))),(RESC)*PARREST.('sigmam')/(150*kk)])
                            cl=F(x0);
                            %if VERBOSE
                                cl
                            %end
                            kk=kk+1;
                        end
                        kk=1;
                        fprintf('Helping p...');
                        while (sum(abs(cl(1:end-1)))>0.4) && (kk<10)
                            A=1/(800);
                            x0=x0 +cl.*(A*[ones(size(x0(1:end-1))),0]);
                            [cl]=F(x0); 
                            %if VERBOSE
                                cl
                            %end
                            kk=kk+1;
                        end  
                        WARNINGS=0;
                        kkk=kkk+1;
                    end
                    %options = optimoptions('fsolve','MaxIter',500,'Display','iter','MaxFunctionEvaluations',5000,...
                    %'FunctionTolerance',10^(0), 'Algorithm', 'trust-region');
                    [output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0,options);
                    %[output,FVAL,EXITFLAG,OUTPUT]= fminsearch(F_ob,x0,options);


                %if (EXITFLAG~=1) && (EXITFLAG~=2)&& (EXITFLAG~=3) && (EXITFLAG~=4)

                    %WARNINGS=WARNINGS+0.5; % had to lower standards
                 %   fprintf('Trying again - with fminsearch');
                    %[output,FVAL,EXITFLAG,OUTPUT]= fsolve(F,x0*0.9,options);
                  %  options = optimset('Display','iter');
                   % options.TolFun=TOL;
                    %options.TolX=TOL;
                    %F_ob=@(x) norm(F(x))^2; 
                    %[output,FVAL,EXITFLAG,OUTPUT]= fminsearch(F_ob,x0,options);
                   

               % end
                if (EXITFLAG~=1) && (EXITFLAG~=2) && (EXITFLAG~=3) && (EXITFLAG~=4)              
                    %WARNINGS=WARNINGS+1;
                    fprintf('Solve failed');
                    FVAL=FVAL
                    OUTPUT
                    EXITFLAG=991;
                 end
                %{
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
                %}
                
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
        