


function [V,workh, workw, Pnn, Pnw0,Pnw, Pw0n,Pw0w0,Pw0w, Pwn, Pww0,Pww,conexp1, conexp2]=...
    locfirst(i,jh,jw,th,tw,lambda,mI1_,mA1_,mIG1_,mAG1_,mI2_,mA2_,mIG2_,mAG2_,vw0w0,vw0n,vnw0,vww,vwn,vnw,vnn,vww0,vw0w)


%todo: check the use of lambda here - somehow high lambda suggests lower
%conexp for husbands???????

%lambda=EQS.('lambda'); 

% barrier to vectorizing? minmax should be ok actually, if I make 0s all a
% matrix
% inputs would be all the v columns as well as stuff to create upper and
% lower bounds of distributions as matrices itself.


%[mI1_,mA1_,mIG1_,mAG1_]=matchdist(i,jh,th,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);
%[mI2_,mA2_,mIG2_,mAG2_]=matchdist(i,jw,tw,mA,mI,mAL,mIL, typeic,D,mm,JLs, betah);


denom=(mA1_-mI1_)*(mA2_-mI2_);
denomG=(mAG1_-mIG1_)*(mAG2_-mIG2_);

if denom*denomG>0
    
%[vbar,vbari]=max([vnw0,vw0n,vw0w0]);
%[vwbar,vwbari]=max([vwn,vww0]);
%[vbarw,vbarwi]=max([vnw,vw0w]);

% not working or taking the global offer if partner accepted local offer
    Pw0n_=max(0,min((vw0n-vw0w)/(1-lambda) - mIG2_,mAG2_-mIG2_))/(mAG2_-mIG2_); % cond on w0n or w0w -about wife
    Pnw0_=max(0,min((vnw0-vww0)/((lambda)) - mIG1_,mAG1_-mIG1_))/(mAG1_-mIG1_); % cond on nw0 or ww0 -about husband
    conexp2__=(1-Pw0n_)*(mAG2_+max((vw0n-vw0w)/((1-lambda)),mIG2_))/2;
    conexp1__=(1-Pnw0_)*(mAG1_+max((vnw0-vww0)/((lambda)),mIG1_))/2;
   
    vwbar=Pw0n_*vw0n + (1-Pw0n_)*vw0w + (1-lambda)*conexp2__;
    vbarw=Pnw0_*vnw0 + (1-Pnw0_)*vww0 + lambda*conexp1__;
   
    db_x=mIG2_*(1-lambda)/lambda - (vwn-vnw)/(lambda); %line intercept
    db_y=mIG1_*(lambda)/(1-lambda) + (vwn-vnw)/((1-lambda));
    db_X=mAG2_*(1-lambda)/lambda - (vwn-vnw)/(lambda);
    db_Y=mAG1_*(lambda)/(1-lambda) + (vwn-vnw)/((1-lambda));
%not working or taking a global offer if partner also has a global offer   
    tresh1=max(min([db_X,(vnw-vww)/(lambda),mAG1_]),mIG1_);
    tresh2=max(min([db_Y,(vwn-vww)/((1-lambda)),mAG2_]),mIG2_);
    tresh1_=min(max([db_x,mIG2_*(1-lambda)/lambda+(vnw-vwn)/((lambda)),mIG1_]),mAG1_); %todo: simplify- this seems like db_X
    tresh2_=min(max([db_y,mIG1_*(lambda)/(1-lambda)+(vwn-vnw)/((1-lambda)),mIG2_]),mAG2_);
    
    MM_=max(0,min(mAG2_-mIG2_,mAG2_-(vwn-vww)/((1-lambda))))...
        *max(0,min(mAG1_-mIG1_,mAG1_-(vnw-vww)/(lambda)))/denomG; 
    Pww_=MM_;
    NN_=max(0,mAG2_-max([mIG2_,tresh2]))*...
        max(0,min(mAG1_-mIG1_,tresh1 - mIG1_ ))/denomG;
    OO_=max(0,min(mAG2_-mIG2_,tresh2 - mIG2_))*...
        max(min(mAG1_-mIG1_,mAG1_-tresh1  ),0)/denomG;
    PP_=0.5*max(0,...
        tresh2-tresh2_)*...
        max(0,tresh1-tresh1_)...
        *(MM_+NN_+OO_<1)/denomG; % does this cover all the zero cases?
    RR_=max(0,...
        min(tresh2_-mIG2_,mAG2_-mIG2_))*...
        max(0,min(tresh1-mIG1_,mAG1_-mIG1_))...
        *(MM_+NN_+OO_<1)*(mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))>=mIG2_)/denomG...
        +...
        max(0,...
        min(tresh2-mIG2_,mAG2_-mIG2_))*...
        max(0,(tresh1_-mIG1_)...
        )*(MM_+NN_+OO_<1)*(mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))<mIG2_)/denomG;
    
    Pwn_=OO_+PP_+RR_*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))>=mIG2_);
    Pnw_=NN_+PP_+RR_*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/((1-lambda))<mIG2_);
    
    %todo: replace tresh1 when appropriate
    conexp1_=(mAG1_+tresh1)*0.5*(OO_+MM_) +...
	     (tresh1+mIG1_)*0.5*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/(1-lambda)>=mIG2_)*RR_+...
         triangle(tresh1_,tresh1,tresh2_,tresh2,0,0)*PP_;
        
    conexp2_=(mAG2_+tresh2)*0.5*(NN_+MM_) +...
	    (tresh2+mIG2_)*0.5*( mIG1_*lambda/(1-lambda)+(vwn-vnw)/(1-lambda)<mIG2_)*RR_+...
        triangle(tresh2_,tresh2,tresh1_,tresh1,0,0)*PP_; %todo - check triangle mAL_kes sense

    vbar=Pww_*(vww) + Pwn_*vwn + Pnw_*vnw + lambda*conexp1_+(1-lambda)*conexp2_;
     
    %todo: probably the above also needs the correction with line
    %intercepts!
    
    dbind=(vbar+vw0w0-vwbar-vbarw>0); 
    db1_x=-mA2_*(1-lambda)/lambda +(vbar-vw0w0)/lambda;
    db1_y=-mA1_*(lambda)/(1-lambda) +(vbar-vw0w0)/(1-lambda);
    db0_x=mI2_*(1-lambda)/lambda - (vwbar-vbarw)/lambda;
    db0_y=mI1_*(lambda)/(1-lambda) + (vwbar-vbarw)/(1-lambda);
    db1_X=-mI2_*(1-lambda)/lambda +(vbar-vw0w0)/lambda;
    db1_Y=-mI1_*(lambda)/(1-lambda) +(vbar-vw0w0)/(1-lambda);
    db0_X=mA2_*(1-lambda)/lambda - (vwbar-vbarw)/lambda;
    db0_Y=mA1_*(lambda)/(1-lambda) + (vwbar-vbarw)/(1-lambda);
    
    db0_tri_x=min(mA1_,max([max((vbar-vwbar)/lambda,db0_x),(vbar-vwbar)/(lambda),mI1_]));
    db0_tri_X=max(mI1_,min([max((vbar-vwbar)/lambda,db0_X),(vbarw-vw0w0)/(lambda),mA1_]));
    db0_tri_y=min(mA2_,max([max((vbar-vbarw)/(1-lambda),db0_y),(vbar-vbarw)/(1-lambda),mI2_]));
    db0_tri_Y=max(mI2_,min([max((vbar-vbarw)/(1-lambda),db0_Y),(vwbar-vw0w0)/(1-lambda),mA2_]));
    
    db1_tri_x=min(mA1_,max([max((vbarw-vw0w0)/lambda,db1_x),(vbarw-vw0w0)/(lambda),mI1_]));
    db1_tri_X=max(mI1_,min([max((vbarw-vw0w0)/lambda,db1_X),(vbar-vwbar)/(lambda),mA1_]));
    db1_tri_y=min(mA2_,max([max(db1_y,(vwbar-vw0w0)/(1-lambda)),(vwbar-vw0w0)/(1-lambda),mI2_]));
    db1_tri_Y=max(mI2_,min([max(db1_Y,(vwbar-vw0w0)/(1-lambda)),(vbar-vbarw)/(1-lambda),mA2_]));
    
    if dbind==1
        MM=max(0,mA2_-max([mI2_,(vwbar-vw0w0)/(1-lambda),db1_y]))*...
            max(0,mA1_-max([mI1_,(vbarw-vw0w0)/(lambda),db1_x]))/denom; 
        % I will be adding a triangle in some cases
        NN=max(0,min(mA2_)-max([mI2_,(vbar-vbarw)/(1-lambda)]))*...
            max(0,min([mA1_,(vbarw-vw0w0)/(lambda)])-max(mI1_))/denom;
        OO=max(0,min([mA2_,(vwbar-vw0w0)/(1-lambda)])-max(mI2_))*...
            max(0,min(mA1_)-max([mI1_,(vbar-vwbar)/(lambda)]))/denom;
        % I will be subtracting a triangle in some cases
        UU=max(0,min([mA2_,(vbar-vbarw)/(1-lambda),db1_Y]) - mI2_)*...
            max(0,min([mA1_,(vbar-vwbar)/(lambda),db1_X])  - mI1_)/denom;

    else
        MM=max(0,mA2_-max(mI2_,(vwbar-vw0w0)/(1-lambda)))*...
            max(0,mA1_-max(mI1_,(vbarw-vw0w0)/(lambda)))/denom; 
        % I will be adding a triangle in some cases
        NN=max(0,min(mA2_)-max([mI2_,db0_y,(vbar-vbarw)/(1-lambda)]))*...
            max(0,min([db0_X,mA1_,(vbarw-vw0w0)/(lambda)])-max(mI1_))/denom;
        OO=max(0,min([db0_Y,mA2_,(vwbar-vw0w0)/(1-lambda)])-max(mI2_))*...
            max(0,min(mA1_)-max([mI1_,(vbar-vwbar)/(lambda),db0_x]))/denom;
        % I will be subtracting a triangle in some cases
        UU=max(0,min(mA2_-mI2_,(vbar-vbarw)/(1-lambda) - mI2_))*...
            max(min(mA1_-mI1_,(vbar-vwbar)/(lambda)  -mI1_),0)/denom;
    end
    SS=0.5*max(0,db1_tri_Y-db1_tri_y)*...
        max(0,db1_tri_X-db1_tri_x)/denom;    % SS ~ dbind==1
    TT=0.5*max(0,db0_tri_Y-db0_tri_y)*...
        max(0,db0_tri_X-db0_tri_x)/denom;
    
    if dbind==1
        E1_1_=max(0,(mA2_-mI2_)*(mA1_-db1_tri_X))/denom; 
        E2_1_=max(0,db1_tri_X-db1_tri_x)*...
            (mA2_-db1_tri_y)/denom;
        E3_=SS;
        conexp1a=(((mA1_+db1_tri_X))/2)*E1_1_+...
            (max(db1_tri_X+db1_tri_x)/2)*E2_1_+...
            triangle(...
            db1_tri_x,db1_tri_X,...
            db1_tri_y,db1_tri_Y,...
            1,1)*E3_; %todo -  check triangle between the two x1 and x2
        E1_2_=max(0,(mA1_-mI1_)*(mA2_-db1_tri_Y))/denom; 
        E2_2_=max(0,db1_tri_Y-db1_tri_y)*...
            (mA1_-db1_tri_x)/denom;
        conexp2a=(((mA2_+db1_tri_Y))/2)*E1_2_+...
            (max(db1_tri_Y+db1_tri_y)/2)*E2_2_+...
            triangle(...
            db1_tri_y,db1_tri_Y,...
            db1_tri_x,db1_tri_X,...
            1,1)*E3_; 
    else
        E1_1_=max(0,(mA2_-mI2_)*(mA1_-db0_tri_X))/denom; 
        E2_1_=max(0,db0_tri_X-db0_tri_x)*...
            (db0_tri_y-mI2_)/denom;
        E3_=TT;
        conexp1a=((mA1_+db0_tri_X)/2)*E1_1_+...
            (max(db0_tri_X+db0_tri_x)/2)*E2_1_+...
            triangle(...
            db0_tri_x,db0_tri_X,...
            db0_tri_y,db0_tri_Y,...
            0,0)*E3_; 
        E1_2_=max(0,(mA1_-mI1_)*(mA2_-db0_tri_Y))/denom; 
        E2_2_=max(0,db0_tri_Y-db0_tri_y)*...
            (db0_tri_x-mI1_)/denom;
        conexp2a=(((mA2_+db0_tri_Y))/2)*E1_2_+...
            (max(db0_tri_Y+db0_tri_y)/2)*E2_2_+...
            triangle(...
            db0_tri_y,db0_tri_Y,...
            db0_tri_x,db0_tri_X,...
            0,0)*E3_; 

    end
    Pw0w0=MM-dbind*SS;
    Pvbar=UU-dbind*SS;
    Pwbar=OO-(1-dbind)*TT;
    Pbarw=NN-(1-dbind)*TT;
    
    V=Pw0w0*(vw0w0) + Pwbar*vwbar + Pbarw*vbarw +Pvbar*vbar+ lambda*conexp1a+(1-lambda)*conexp2a;

    conexp1=conexp1a+Pvbar*conexp1_+Pbarw*conexp1__;
    conexp2=conexp2a+Pvbar*conexp2_+Pwbar*conexp2__;
    %should the other stuff have probabilities in front of them?
 
    Pnw=Pvbar*Pnw_;
    Pwn=Pvbar*Pwn_;
    Pww=Pvbar*Pww_;

    Pw0n=Pwbar*Pw0n_;
    Pw0w=Pwbar*(1-Pw0n_);

    Pnw0=Pbarw*Pnw0_;
    Pww0=Pbarw*(1-Pnw0_);
    
    Pnn=0;
    workh=Pww+Pwn+Pw0w0+Pw0n+Pw0w+Pww0;
    workw=Pww+Pnw+Pw0w0+Pnw0+Pw0w+Pww0;  
    
    
    
check=round(Pww+Pwn+Pw0w0+Pw0n+Pw0w+Pww0+Pnw+Pnw0,2);
if (round(Pw0w0+ Pwbar+ Pbarw+Pvbar,2)==1)*(round(Pww_ + Pwn_+ Pnw_,2)==1)...
*(check==1)*(Pww<=1)*(Pvbar<1 | (conexp1a==0 & conexp2a==0))...
~=1
    check
    (Pw0w0+ Pwbar+ Pbarw+Pvbar)
    (Pww_ + Pwn_+ Pnw_)
    Pww<=1
    %(conexp1<max(mA1_,mAG1_))*(conexp1>=min(mI1_,mIG1_)*workh)*(conexp2<max(mA2_,mAG2_))*(conexp2>=min(mI2_,mIG2_)*workw)
    %(conexp1<max(mA1_,mAG1_))
    %(conexp1>=min(mI1_,mIG1_)*workh)
    %(conexp2<max(mA2_,mAG2_)) % I think this does not work with the match shocks being negative on average?
    %(conexp2>=min(mI2_,mIG2_)*workw)
    (Pvbar<1 | (conexp1a==0 & conexp2a==0))
    dbind
    i
    jh
    jw
    fprintf('Warning')
    
    % I seem to have issues at very low disperions
    return
end

else % this is not really correct, because I do not have the conexp here.... I only use it to judge speed of above
    Pwn=0;
    Pww=0;
    Pnw=0;
    Pww0=0;
    Pw0w=0;
    if vww>vwn && vww>vnw
        Pww=1;
        V=vww;
        workh=1;
        workw=1;
    elseif vwn>vww && vwn>vnw
        Pwn=1;
        V=vwn;
        workh=1;
        workw=0;
    elseif vnw>vww && vnw>vwn
        Pnw=1;
        V=vnw;
        workh=0;
        workw=1;
    elseif vnw==vww && vnw>vwn
        Pnw=0.5;
        Pww=0.5;
        V=vnw*Pnw+vww*Pww; 
        workh=0.5;
        workw=1;
    elseif vwn==vww && vwn>vnw
        Pwn=0.5;
        Pww=0.5;
        V=vwn*Pwn+vww*Pww; 
        workh=1;
        workw=0.5;
    elseif vnw==vwn && vnw>vww
        Pnw=0.5;
        Pwn=0.5;
        V=vnw*Pnw+vwn*Pwn; 
        workh=0.5;
        workw=0.5;
    end
   conexp1=0;
   conexp2=0;
   Pnw0=0;
   Pw0n=0;
   Pw0w0=0;
   Pnn=0;
end




end