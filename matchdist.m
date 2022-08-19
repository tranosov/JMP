function [mIL1_,mAL1_,mI1_,mA1_,L_low,low]=matchdist(i,j,t,mA,mI,mAL,mIL,type,D,mm,JLs, betah) % local first!
global Davg ICLOWS
if isempty(mm)
    mm=0;
end
type=round(type,10);

% just one suburb bad:
if type==1
    lowest_t1=find(JLs(:,1)==min(JLs(:,1)));
    lowest_t2=find(JLs(:,2)==min(JLs(:,2)));
    % local: 
    L_low =((i==lowest_t1 & t==1) | (i==lowest_t2 & t==2)); % do I live in my hub?
    % sector 2 is female, neighb 3 (first suburb) is dominated
    low   =((j==lowest_t1 & t==1) | (j==lowest_t2 & t==2));

    %L_low =((i~=1)); % do I live in my hub?
    %low   =((j~=1));

    if L_low
        mIL1_=mIL;
        mAL1_=mAL;
    else
        mIL1_=mI;
        mAL1_=mA;
    end


    if low
        mI1_=mIL;
        mA1_=mAL;
    else
        mI1_=mI;
        mA1_=mA;
    end
end
if type==2 % slide
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    low=betah*d(j)+1-betah*10;
    
    mIL1_=mI;
    mAL1_=mA;
    mI1_=mI;
    mA1_=mA;
    L_low=betah*d(i)+1-betah*10;
    
end

if type==3 % slide + for the constant as well!
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    low=betah*d(j)+1-betah*10;
    
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low=betah*d(i)+1-betah*10;
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end

if type==3.5 % slide + for the constant as well - to be above 0!
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    low=betah*d(j)+1-betah*19.1400;
    
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low=betah*d(i)+1-betah*19.1400; % make it so it is all positive here... but it totally arbitrary. longer distances - all matches worse!
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end

if type==3.6 % slide + max
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    low=betah*d(j)+1-betah*max(d); % characteristic - longer distance - city has BETTER matches, suburbs same
    
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low=betah*d(i)+1-betah*max(d); % make it so it is all positive here...
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end


if (type==3.7) || (type==7) % slide + average  - for hours can also be filled in with an exp!
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    davg= sum(JLs(:,t).*d');
    low=betah*d(j)+1-betah*davg; % average 'productivity' stay the same. but issue - is negative ic in suburbs. bad for husbands
    
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low=betah*d(i)+1-betah*davg; % make it so it is all positive here...
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end

if (type==3.71) || (type==71) % slide + average  - for hours can also be filled in with an exp!
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    if isempty(Davg) 
        davg= sum(JLs(:,t).*d');
        Davg=zeros(2,1);
        Davg(t)=davg;
    elseif Davg(t)==0
        davg= sum(JLs(:,t).*d');
        Davg(t)=davg;
    else
        davg=Davg(t);
    end
    davg= sum(JLs(:,t).*d');
    low=betah*d(j)+1-betah*davg; % average 'productivity' stay the same. but issue - is negative ic in suburbs. bad for husbands
    
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low=betah*d(i)+1-betah*davg; % make it so it is all positive here...
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end

if type==6 % slide + average + super positive!
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    davg= sum(JLs(:,t).*d');
    low=betah*d(j)-betah*davg; % average 'productivity' stay the same. 
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low=betah*d(i)-betah*davg; % make it so it is all positive here...
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end


if type==4
    const=500;
    I=size(JLs);
    I=I(1);
    %T=s(2);
    for u=1:I
        d(u)=Do(u,JLs(:,t),D);
    end
    low= (1-exp(const*(-betah*d(j)))+exp(const*(-betah*10))); %betah*d(j)+1-betah*10;
    
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    L_low= (1-exp(const*(-betah*d(i)))+exp(const*(-betah*10))); %betah*d(i)+1-betah*10;
    
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
    if mm~=0
        fprintf('I think mm should be 0.')
    end
   
end


if (type==99) % FIXED!
    %I=size(JLs);
    %I=I(1);
    %T=s(2);
    %for u=1:I
    %    d(u)=Do(u,JLs(:,t),D);
    %end
    %davg= sum(JLs(:,t).*d');
    %low=betah*d(j)+1-betah*davg; % average 'productivity' stay the same. but issue - is negative ic in suburbs. bad for husbands
    low=ICLOWS(j,t);
    
    mI1_=mI+mm*(1-low);
    mA1_=mA+mm*(1-low);
    %L_low=betah*d(i)+1-betah*davg; % make it so it is all positive here...
    L_low=ICLOWS(i,t);
    mIL1_=mI+mm*(1-L_low); %+ic*mm
    mAL1_=mA+mm*(1-L_low);
end

%{
mIL1_=0;
mAL1_=0;
mI1_=0;
mA1_=0;
%}
end

%todo: so that low is not 1-0 but is 0 for center, between 1 and 0 for
%suburb hub and 1 for wrong suburb. make the scaler be the distance to an
%average job - if distances increase, things that are not a hub are even
%worse!