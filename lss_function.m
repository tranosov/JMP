
function [lssh,lssw,p,LA,time,EXITFLAG]...
    =lss_function(MtoF,params,EQS,PARREST,psolve,msolve)

global VERBOSE

    LA0=params{'LA0',:};
    p0=[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:}];
    params('MtoF','value')={MtoF};
    reinitialize=0.5; % 0 vs 0.5 vs 1 - 0.5 is reset PARREST and EQS with new lambda, but not new inputs  
    [p,LA,EXITFLAG,time]=solvemodel(params,EQS,PARREST,p0,LA0,reinitialize,psolve,msolve);
    if EXITFLAG==999
        lssh=999;
        lssw=999;
        p=p0;
        LA=LA0; 
        return
    else
        tic
        HS=PARREST.('HS');
        Jw=PARREST.('Jw');
        sigmam=PARREST.('sigmam');
        I=size(HS,2);
        T=size(Jw,2);
        W=3;
        
        
        [DS,~, DSh, DSw, VS,~]=DSingle(p,0,EQS,PARREST);
        if (size(DS,1)==1)
            lssh=999;
            lssw=999;
            EXITFLAG=999;
            return
        end
        [DC, DC1, DC2,~,Pw,~,~,cexp1,cexp2,OUTC]=DCouple(p,1,LA,EQS,PARREST); 
        if (size(DC,1)==1)
            lssh=999;
            lssw=999;
            EXITFLAG=999;
            return
        end

        DSh_agr=sum(DSh,4);
        DSw_agr=sum(DSw,4);
        uhS=2*sum(sum(sum(DSh_agr.*VS)))./sum(sum(sum(DSh_agr)));
        uwS=2*sum(sum(sum(DSw_agr.*VS)))./sum(sum(sum(DSw_agr)));

        uh=OUTC.('uh') ;
        uw=OUTC.('uw') ;
        uhC_det=zeros(I,I,T,T,I,W,W); 
        uwC_det=zeros(I,I,T,T,I,W,W);
        uhC=zeros(I,I,T,T,I); 
        uwC=zeros(I,I,T,T,I);
        for jh=1:I
            for jw=1:I
                for th=1:T
                    for tw=1:T
                        for i=1:I
                            for wh=1:3 
                               for ww=1:3
                                   hi=wh>1;
                                   wi=ww>1;
                                   wwi=hi*(1-wi)+2*wi*(1-hi) + 3*wi*hi;
                                   jh_=jh*(wh==3) + i*(wh==2|wh==1);
                                   jw_=jw*(ww==3) + i*(ww==2|ww==1);
                                   if wwi>0
                                   uhC_det(jh,jw,th,tw,i,wh,ww)=uh(th,tw,jh_,jw_,i,wwi) ; 
                                   uwC_det(jh,jw,th,tw,i,wh,ww)=uw(th,tw,jh_,jw_,i,wwi) ;
                                   end
                                end
                            end      
        uhC(jh,jw,th,tw,i)= sum(sum(Pw(jh,jw,th,tw,i,:,:).*uhC_det(jh,jw,th,tw,i,:,:))) + cexp1(jh,jw,th,tw,i);
        uwC(jh,jw,th,tw,i)= sum(sum(Pw(jh,jw,th,tw,i,:,:).*uwC_det(jh,jw,th,tw,i,:,:))) + cexp2(jh,jw,th,tw,i);

                        end
                    end
                end
            end
        end

        DC1_agr=sum(sum(DC1,7),6);
        DC2_agr=sum(sum(DC2,7),6);

        uhC=sum(sum(sum(sum(sum(DC1_agr.*uhC)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uhC)))))./sum(sum(sum(sum(sum(DC2_agr)))));
        uwC=sum(sum(sum(sum(sum(DC1_agr.*uwC)))))./sum(sum(sum(sum(sum(DC1_agr)))))+sum(sum(sum(sum(sum(DC2_agr.*uwC)))))./sum(sum(sum(sum(sum(DC2_agr)))));
        ssh=1/(1+exp((uhC-uhS)/sigmam));
        ssw=1/(1+exp((uwC-uwS)/sigmam)); 
        lssh=log( (1/ssh) -1);
        lssw=log( (1/ssw) -1);

        time=time+toc;
        if VERBOSE
            toc
        end
    end
end