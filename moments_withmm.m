
function [moments_,other,time,EXITFLAG,params]...
    =moments_withmm(p,LA,params,EQS,PARREST,withmm,fast,FORCEFIT)
global RESC VERBOSE
tic
[DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC, DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS,EQS,PARREST] = ...
    populations(p,LA,EQS,PARREST,FORCEFIT);
other=99;
[moments_,~]...
                =moments(p,EQS,PARREST, DC, DC1, DC2,HCouple,Pw,hC,VC,cexp1,cexp2,OUTC, DS,HSingle, DSh, DSw, VS,Pws,cexps,OUTS,fast);
time=toc;
if VERBOSE
toc
end


if withmm
    %update p0, LA0
    EQS.('inputs')=OUTC.('inputs');
    EQS.('inputsS')=OUTS.('inputs'); % to be sure stuff gets initialized with the  last solution

    clear global IN INS % why?

    params=PARREST.('params'); % has the updated thetas!- and also to output + updated prices
    params('LA0','value')={LA};
    params('p0_1','value')={p(1)};
    params('p0_2','value')={p(2)};
    params('p0_3','value')={p(3)};
    %params('sstaysingleh','value')={LA};
    %params('sstaysinglew','value')={LA};
    mf1=exp(0.01);
    mf2=exp(-0.01);
    %{
    [lssh,lssw,p0,LA0]...
    =lss(1,params,0,1)
    [lssh1,lssw1,p1,LA1]...
    =lss(mf1,params,0,1)
    [lssh2,lssw2,p2,LA2]...
    =lss(mf2,params,0,1) % down - seems to be a more drastic change?
    L=1/(1-(lssh1-lssh)/(lssw1-lssw))
    L=1/(1-(lssh-lssh2)/(lssw-lssw2))
    L=1/(1-(lssh1-lssh2)/(lssw1-lssw2)) % somehow this is closer to the truth if  the dif is a bit bigger?

    derw=(lssw1-lssw2)/(mf1-mf2) % way too much - not sure why?
    derh=(lssh1-lssh2)/(mf1-mf2)
    %}
    %F=@(x) Clearing_withmm(x,EQS,PARREST); 
    %F([log(p),LA*(RESC)])
    
    global mup
    mup=1;
    [lssh1_,lssw1_,p1_,LA1_,time_,EXITFLAG]...
    =lss_function(mf1,params,EQS,PARREST,1,1);
    time=time+time_;
    if EXITFLAG==999
        return
    end
    
    mup=0;
    [lssh2_,lssw2_,p2_,LA2_,time_,EXITFLAG]...
    =lss_function(mf2,params,EQS,PARREST,1,1);
    time=time+time_;
    if EXITFLAG==999
        return
    end
        %[lssh,lssw,p0_,LA0_]...
        %=lss(1,params,1,1);


    %L=1/(1-(lssh1_-lssh)/(lssw1_-lssw))
    %L=1/(1-(lssh-lssh2_)/(lssw-lssw2_))
    L=1/(1-(lssh1_-lssh2_)/(lssw1_-lssw2_)); % so far - seems almost exactly as above, which is good

    %derw=(lssw1_-lssw2_)/(mf1-mf2) % WAAAAAY too big. like omg
    %derh=(lssh1_-lssh2_)/(mf1-mf2)
    %(LA1_+LA2_)/2

moments_('L',:)={L};
end

end
