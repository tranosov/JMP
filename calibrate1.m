

function parsc=calibrate1(params,momentall)

b=params{'b',:};
%ce_=params{'ce_',:};
%crra_=params{'crra_',:};
l=params{'LA0',:};
wa=params{'wa',:};
wgap_raw=params{'wgap_raw',:};
crrat_=params{'crrat_',:};
wc=params{'wc',:};
pid_=params{'pid_',:};

ic0_=params{'ic0_',:};
Xi=params{'Xi',:};
crrah_=params{'crrah_',:};
crra_=params{'crra_',:};

typeic=params{'typeic',:};
if typeic==7
    fxi=@(a) exp(a);
else
    fxi=@(a) a;
end

a=1;
ss_=momentall{'snmarried',:};
ss=(ss_*(2/3)+(1/3));
spww=1+(momentall{'wlfp_dif',:});

ls0=momentall{'shours',:}/(24*365);
d0=momentall{'scommiles',:};
xs0=momentall{'shwk',:}/(24*365); % housework of singles
%Ls0=1-a*ls0 -b*d0-xs0;

lsh0=momentall{'hhours_pww',:}/(24*365);
lsw0=lsh0+momentall{'whours_pww_dif',:}/(24*365);

xw0=momentall{'whwk_pww',:}/(24*365);
xh0=xw0+momentall{'hhwk_pww_dif',:}/(24*365);
dh0=momentall{'scommiles',:} - momentall{'shcommiles_dif',:};
dh0_h=dh0+0.5; %this is a bit arbitrary

dw0=momentall{'scommiles',:} - momentall{'swcommiles_difw',:}; % this is not totally fair - probably the gap is smaller among both work. couples
lsh0_h=lsh0+momentall{'hhours_pwn_dif',:}/(24*365);
xh0_h=xh0 + momentall{'hhwk_pwn_difhpww',:}/(24*365); % corrected
xw0_h=xw0+momentall{'whwk_pwn_dif',:}/(24*365);
Lh0=1-lsh0-xh0-b*dh0;
Lh0_h=1-lsh0_h-xh0_h-b*dh0_h;
Lw0=1-lsw0-xw0-b*dw0;
Lw0_h=1-xw0_h;
w1=wa*exp(ic0_*wc);
w2=wa*exp(wgap_raw)*exp(ic0_*wc);

%crrat_=(log(w1/w2)+log((1-l)/l))/log((1-lsw0-xw0-b*dw0)/(1-lsh0-xh0-b*dh0) ) ; % take out of
%the calibarion - it is too important!
piel_=crrat_*log(((Lw0_h)*(Lh0))/((Lw0)*(Lh0_h)))/log((xh0*xw0_h)/(xh0_h*xw0)) ;
pih_=((l/(1-l))*(xh0/xw0)^piel_)/((l/(1-l))*(xh0/xw0)^piel_+((Lh0)/(Lw0))^crrat_ );
piw_=1-pih_;
T0= (piw_.*xw0.^(1-piel_)+ (1-piw_).*xh0.^(1-piel_) ).^(1/(1-piel_));
T0_h=(piw_.*xw0_h.^(1-piel_)+ (1-piw_).*xh0_h.^(1-piel_) ).^(1/(1-piel_));

crrax_=(crrat_*log((Lh0_h)/(Lh0))+piel_*log(T0_h/T0)-piel_*log(xh0_h/xh0) )/(log(T0_h/T0) );
pi_= (l/pih_).*(T0^(crrax_-piel_))*(xh0^(piel_))/((Lh0)^crrat_); %TR: corrected! + treat as relative to value of time
pish_=(xs0^(crrax_))/((1-a*ls0 -b*d0-xs0)^crrat_);

NLY_=momentall{'NLY_',:};
NLY=( (1-ss)*((lsh0*w1+lsw0*w2)*spww  + (1-spww)*(lsh0_h*w1)) + ss*ls0*w1*2   )*NLY_/(1-NLY_);
% depends on wage gap. and on wc

%(lsh0*w1+lsw0*w2)*NLY_/(1-NLY_);

piw_=piw_-b*(dh0-dw0)*pid_;




w1_exp=@(ls)wa.*exp(wc*ic0_);
wd1=@(ls)wa.*exp(wc*ic0_);
w2_exp=@(ls)wa.*exp(wgap_raw+wc*ic0_);
Y0=(lsh0*w1_exp(lsh0)+lsw0*w2_exp(lsw0)) + NLY;
mucouple=((l^(1/crra_)+(1-l)^(1/crra_))^crra_)/((Y0- 2)^crra_)
le_=(1/l)*(mucouple*(wd1(lsh0))+fxi(Xi*ic0_))*(1-a*lsh0 -b*dh0-xh0)^crrat_


mucouple_check=(1/(wd1(lsh0)))*((l*le_)/(1-a*lsh0 -b*dh0-xh0)^crrat_-fxi(Xi*ic0_)) 
if mucouple<0
    fprintf('IMPLIES IMPOSSIBLE MUS')
end

eta_=2*mucouple*(1)^(crrah_);
%ceh_=(mucouple*(Y0- 2)^crra_)/((l^(1/crra_)+(1-l)^(1/crra_))^crra_) ;



parsc = array2table([NLY,pish_,piel_,crrax_,piw_,pi_]','VariableNames',{'value'}, 'RowNames',...
    {'NLY','pish_','piel_','crrax_','piw_','pi_'}); %'crrat_',log(ce_) %'ces_','crrah_'

end