%{

%}

function [fid]=Table_main(moments_,momentest)

%% Load data from Table

%cd '/Users/jablanco/Dropbox/Papers/Nominal Anchor/data/output/Tables'
cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'

%data_m  =xlsread('C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Output\tables\Moments.xlsx',1,'B1:B21');
%Data_load_m  =[Data_load_m(:,1:4),Data_load_m(:,9:end),Data_load_m(:,5:8)];


Keys = {'Average commute of single $d^s$','$d^s_h-d^h$','$d^s_w-d^w$',...
    '$h^w_{\\textrm{both work}}-h^h_{\\textrm{both work}}$','$h^h_{\\textrm{both work}}$', '$h^h_{\\textrm{just husband works}}-h^h_{\\textrm{both work}}$', '$h^w_{\\textrm{just wife works}}-h^w_{\\textrm{both work}}$',...
    '$h^s$','$x^w_{\\textrm{both work}}$', '$x^w_{\\textrm{just husband works}}-x^w_{\\textrm{both work}}$','$x^w_{\\textrm{just wife works}}-x^w_{\\textrm{both work}}$',...
    '$x^h_{\\textrm{both work}}-x^w_{\\textrm{both work}}$','$x^h_{\\textrm{just husband works}}-x^h_{\\textrm{both work}}$','$x^h_{\\textrm{just wife works}}-x^h_{\\textrm{both work}}$',...
    '$x^s$', 'LFP of wives-husbands', 'LFP of husbands', '$\\log(\\frac{w^w}{w^h}) $',...
    'Distance to an average job for a couple ($d_j^h$)','Distance to an average job in own labor market for a husband ($d_o^h$)',...
    'Distance between 2 random jobs','Distance between 2 random jobs of the husbands labor market',...
    'Distance between 2 random jobs, in h and w labor markets', 'Distance between husbands and wifes actual jobs',...
    '$P(city|couple)-P(city|single)$',...
    '$d_o^s - d_o^h$','$d_o^w - d_o^h$','$Var(d_o)$','$\\lvert d_o^w - d_o^h \\lvert$', ...
    '$\\beta_{wd}^{comm}$','$\\beta_{wd}^{work}$','$\\beta_{wd}^{hours+}$','$\\beta_{wd}^{hours}$','$\\beta_{wd}^{x}$',...
    '$\\beta_{d}^{comm}$','$\\beta_{d}^{work}$','$\\beta_{d}^{hours+}$','$\\beta_{d}^{hours}$','$\\beta_{d}^{x}$',...
    '$\\beta_{wd}^{\\log(w)}$',...
    'Share of non-labor income','Expenditure share on housing','Share of population in city','Share of jobs in city',...
    '$\\log(p)$ distance to jobs gradient','$\\log(p)$ city over suburb',...
    '$\\hat{\\lambda}_0$','Share never married'};
    
%additions = {' ',' ',' ',' ',' annual hours',' miles',' miles',' miles',' miles','',' annual hours','','',' \\%',' \\%','','','','','',''};
%prefs = {' ',' ',' ','\~ ',' ',' ',' ',' ',' ','','\~','','',' ',' ','','','','','',''};

moments_.Keys = Keys';
moments_.data = momentest.value;
moments_ = sortrows(moments_,'sort');

Keys=moments_.Keys;
values=table2array(moments_(:,'moments_'));
data=table2array(moments_(:,'data'));
calibrations=table2array(moments_(:,'directlyused'));
%% Print table in Latex


cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Output\tables\'

fid=fopen('MomentsEstimation.tex','w+');
% Declare variables
fprintf(fid,'\\begin{longtable}[htp]{ c | c | c | c } \n');
%fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\footnotesize \n');
%fprintf(fid,'\\caption{Moments \\label{main:tab:businessCycleStatMain}} \n');
%fprintf(fid,'\\begin{tabular}{ c | c | c | c } \n');
%fprintf(fid,' \\\\ \n');
%fprintf(fid,' & \\multicolumn{4}{c|}{} & \\multicolumn{3}{c|}{Correlation of} &  \\multicolumn{4}{c}{ }  \\\\   \n');
fprintf(fid,' Moment & Model value & Data value  & directly used in calibration\\\\ \\cline{1-4}\\cline{1-4}  \n');
%fprintf(fid,' & Long-term &           &        & Nominal       &           &        & Nominal       & Long-term &           &        & Nominal \\\\    \n');
%fprintf(fid,' & Inflation & Inflation & Output & exchange rate & Inflation & Output & exchange rate & Inflation & Inflation & Output & exchange rate \\\\ \\hline \\hline  \n');

%fprintf(fid,'& &&&& &&& &&&  \\\\  \n');


for jj=1:size(Keys,1)
    
     fprintf(fid,Keys{jj,1});
     fprintf(fid,' ');
     fprintf(fid,strcat(' & ',num2str(values(jj,1),'%5.3f')));
     fprintf(fid,' ');
      fprintf(fid,strcat(' & ',num2str(data(jj,1),'%5.3f')));
     fprintf(fid,' ');
     fprintf(fid,strcat(' & ',num2str(calibrations(jj,1),'%5.0f')));
     fprintf(fid,' ');
     %fprintf(fid,'& ')
     %fprintf(fid, additions{1,jj});

     fprintf(fid,'\\\\  \n');
    
    
end
 
%fprintf(fid,'& &&&& &&& &&&  \\\\  \\hline \\hline\n');
%fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\vspace{0.1in} \\\\  \n');
%fprintf(fid,'\\end{center} \n');
%fprintf(fid,'\\justify \n');
fprintf(fid,'\\caption{ \\footnotesize  ');
fprintf(fid,'Moments data versus model - used in estimation. In addition, I use a ratio of average commuting time and distance in miles as a scaling factor, I constarint housing prices to be one on average, and I impose that the ratio of men and women in the metro-area is equal to one.} \\label{tab:moms} ');
 fprintf(fid,'  \n');
fprintf(fid,'\\end{longtable} \n');
fclose(fid);

% todo: figure out where to put untargeted moments? the initial 'gender
% gaps & distance from center moments - present??' and work on enclave to
% make them closer to the model setup.

cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab\'


%movefile('MomentsMain.tex','C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Output\tables\');
end  
