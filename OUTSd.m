function OUTS=OUTSd(OUTS,EQS,PARREST)
%global EQS JLs D AS typeic WARNINGS mm params


typeic=PARREST.('typeic');
D=PARREST.('D');
mm=PARREST.('mm');
JLs=PARREST.('JLs');
AS=PARREST.('AS');
betah=PARREST.('betah');
params=PARREST.('params');
pp=[params{'p0_1',:},params{'p0_2',:},params{'p0_3',:} ]; 
N=size(JLs);
I=N(1);
T=N(2);
wfh=PARREST.('wfh');
T0=T;
if wfh>0
    T=T0*2;
end

us=EQS.('us');
utills=EQS.('utills') ;
utilcls=EQS.('utilcls') ;
Ys=EQS.('Ys');
utilcs=EQS.('utilcs');

vs=zeros(T,I,I,2);
uls=zeros(T,I,I);
ucls=zeros(T,I,I);
ics=zeros(T,I,I);
Yss=zeros(T,I,I);
ucs=zeros(T,I,I);
conss=zeros(T,I,I);


conss_=OUTS.('conss');
Hs=OUTS.('Hs');
lss=OUTS.('lss');
xs=OUTS.('xs');
Yss_=OUTS.('Ys');

for j=1:I
    for i=1:I
        for t=1:T
            d=D(i,j);
            p=pp(i);
            if t>T0
                d=0;
                [~,~,~,~,~,low]=matchdist(i,j,t-T0,0,0,0,0,typeic,D,mm,JLs,betah);
            else
                [~,~,~,~,~,low]=matchdist(i,j,t,0,0,0,0,typeic,D,mm,JLs,betah);
            end
            ic=1-low;
            ics(t,j,i)=ic;
            
            
            conss(t,j,i)=conss_(t,j,i)+(-Yss_(t,j,i)+Ys(lss(t,j,i),ic)); % if there is loss in income - put it on consumption
            Yss(t,j,i)=Ys(lss(t,j,i),ic);
            vs(t,j,i,2)= us(1,d,AS(i),conss(t,j,i),Hs(t,j,i),lss(t,j,i),xs(t,j,i),0,ic);
            uls(t,j,i)=utills(lss(t,j,i),xs(t,j,i),d);
            ucls(t,j,i)=utilcls(d);
            ucs(t,j,i)=utilcs(conss(t,j,i));
            
            

        end
    end
end


OUTS.('vs') = vs;
OUTS.('ics') = ics;
OUTS.('uls') = uls;
OUTS.('ucls') = ucls;
OUTS.('ucs') = ucs;
OUTS.('conss') = conss;
OUTS.('Ys') = Yss;
end