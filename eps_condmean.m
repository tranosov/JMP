function [OUTC,OUTS]=eps_condmean(OUTC,OUTS,PARREST,VS,VC)
sigmal=PARREST.('sigmal');
I=size(VS,1);
T=size(VS,2);
PS_=exp(VS./repmat(sigmal,I,T,I))./repmat(sum(exp(VS./repmat(sigmal,I,T,I)),3),1,1,I);
PC_=exp(VC./(repmat(sigmal*2,I,I,T,T,I)))./repmat(sum(exp(VC./(repmat(sigmal*2,I,I,T,T,I))),5),1,1,1,1,I);

S=10000;
rng(357)
dist=makedist('gev',0,sigmal,0); % gev not ev!
epses = random(dist,[S,I]);

epsS=zeros(I,T,I);
for j=1:I
    for t=1:T
        [values,is]=max(reshape(VS(j,t,:),1,I)+epses,[],2);
        %[sum(is==1)/S,sum(is==2)/S,sum(is==3)/S] should be PS_(j,t,:)
        % so far I have some issue here?
        for i=1:I
            epsS(j,t,i)=sum(epses(is==i,i))/sum(is==i);           
        end
    end
end

epsC=zeros(I,I,T,T,I);
for jh=1:I
    for jw=1:I
        for th=1:T
            for tw=1:T
                [values,is]=max(reshape(VC(jh,jw,th,tw,:),1,I)+2*epses,[],2);
                for i=1:I
                    epsC(jh,jw,th,tw,i)=sum(epses(is==i,i))/sum(is==i);  
                end
            end
        end
    end
end                    

OUTC.('epsC') = epsC;
OUTS.('epsS') = epsS;
end