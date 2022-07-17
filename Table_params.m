%{
"Put title here", put names

Note: This code generates the table II for the paper: "Business Cycle Moments: Model and Data"

%}

function [fid]=Table_params(params_)

%% Load data from Table

%cd '/Users/jablanco/Dropbox/Papers/Nominal Anchor/data/output/Tables'
cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab'

%data_m  =xlsread('C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Output\tables\M.xlsx',1,'B1:B21');
%Data_load_m  =[Data_load_m(:,1:4),Data_load_m(:,9:end),Data_load_m(:,5:8)];



Keys = { '$w_a$','$b$','$\\Omega_c $','$Y$','$\\frac{\\Omega_H}{\\Omega_c}$','$\\omega_c$', '$\\omega_l$','$\\Omega_x$',...
    '$\\bar{\\kappa}_w$', '$\\eta_x$','$\\omega_x$', '$\\Omega_x^s$','$w_{gap}$',...
    '$\\kappa_d$','$\\phi$', '$\\rho$', '$\\delta_w$',...
    '$D(1,2)=D(1,3)$','$D(3,2)/D(1,2)$',...
    '$A(1)$','$A^c(3)$','$\\sigma_{\\epsilon_i}$','$A^c(2)-A^c(3)$',...
    '$f(2|1)/f(3|1)=f(3|2)/f(2|3)$','$f(1)/f(2)$',...
    '$\\pi$','$\\frac{\\bar{N}^C_{1,1}+\\bar{N}^C_{1,2}}{\\sum_{u,v} \\bar{N}^C_{u,v}}$',...
    '$E(\\xi_0)$','$Var(\\xi_0)$','$\\bar{\\Xi}$','$w_{\\Xi}$',...
    '$\\lambda$', '$\\sigma_{\\theta_i}$','$\\Theta_w$','$\\Theta_h-\\Theta_w$','$\\bar{h}^U$'}
    
%{
'$\\beta $','$\\omega_c $','$\\omega_l $','$\\omega_H$','$\\omega_x$','$\\Omega_c$','$\\Omega_H$','$\\Omega_x^s$','$\\Omega_x^c$',...
    '$\\kappa_w^0$','$\\kappa_d$','$\\eta_x$',  '$\\phi$', '$\\rho$', '$\\delta_w$',...
    '$a$','$b$','$\\bar{\\Xi}$','$E(\\xi_0)$','$Var(\\xi_0)$','$\\pi$', '$f(1|1)=f(1|2)$','$f(2|1)/f(3|1)=f(3|2)/f(2|3)$','$\\frac{\\bar{N}^C_{1,1}+\\bar{N}^C_{1,2}}{\\sum_{u,v} \\bar{N}^C_{u,v}}$',...
    '$D(1,2)=D(1,3)$','$D(3,2)/D(1,2)$',...
    '$\\sigma_{\\epsilon_i}$','$A^s(1)$','$A^c(2)=A^c(3)$',...
    '$\\lambda$','$N^S/N in period 2,3$'};
%}    

params_.Keys = Keys';
params_ = sortrows(params_,'sort');

Keys=params_.Keys;
values=table2array(params_(:,'value'));
calibrations=table2array(params_(:,'tocalibrate'));
%{
additions = {'Scaling of commuting distance to time',...
             'To match small hours increases with partner lfp=0',...
             'Hours gaps and wage gap of $5$ percent in couples',...
             'Reasonably inelastic housing demand',...
             'Housework times when both work versus just husband',...
             'Price and quantity of housing 1 on avg', ... % I need it, but not sure why
             ' $35$ percent of income spent on hoursing',...
             'Housework hours by singles',...
             'Housework hours by couples',...
             'Housework hours by men vs women in couples',...
             'Hours reacting to commuting',...%sharpens
             'Leisure and housework of men and women, when both work vs just husband',...
             'Commute of singles',...
             'Commute of singles vs couples',...
             'Commute of men vs women in couples',...
             'Singles hours',... % normalization
             'Wage gap in couples of $5$ percent ' ,...
             'hours of wives more sensitive to being far away from industry hub',...
             'Both men and women inframarginal on participation + LFP of women',...
             '',...
             'Commute of singles (2)',... % this is a second moment I am saying this, which seems to be an issue
             'Share of jobs less than 10 miles away from center',...
             '',...
             '',...
             'Distance to an average job',...
             'Distance to an average job for subgroups',...
             'Distance to a job in own industry',...
             'Singles $7$ percent more likely to live less than 10 miles away from center',...
             '',...
             'Distance to a job in own industry for men vs women in couples',...
             'appr. share of never married adults in metro areas'};
 
%}

%% Print table in Latex


cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Output\tables\'

fid=fopen('ParamsMain.tex','w+');
% Declare variables
fprintf(fid,'\\begin{table}[htp]\n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\footnotesize \n');
%fprintf(fid,'\\caption{Moments \\label{tab:params} \n');
fprintf(fid,'\\begin{tabular}{ c | c | c  } \n');
%fprintf(fid,' \\\\ \n');
%fprintf(fid,' & \\multicolumn{4}{c|}{} & \\multicolumn{3}{c|}{Correlation of} &  \\multicolumn{4}{c}{ }  \\\\   \n');
fprintf(fid,' Parameter &  value & calibrated  \\\\ \\cline{1-3}\\cline{1-3}  \n');
%fprintf(fid,' & Long-term &           &        & Nominal       &           &        & Nominal       & Long-term &           &        & Nominal \\\\    \n');
%fprintf(fid,' & Inflation & Inflation & Output & exchange rate & Inflation & Output & exchange rate & Inflation & Inflation & Output & exchange rate \\\\ \\hline \\hline  \n');

%fprintf(fid,'& &&&& &&& &&&  \\\\  \n');


for jj=1:size(Keys,1)-4
    
     fprintf(fid,Keys{jj,1});
     fprintf(fid,' ');
     fprintf(fid,strcat('& ',num2str(values(jj,1),'%5.3f')));
     
     %fprintf(fid,strcat('& ',num2str(data_m(jj,1),'%5.3f')));
     fprintf(fid,' ');
     fprintf(fid,strcat('& ',num2str(calibrations(jj,1),'%5.0f')));
     fprintf(fid,' ');
     %fprintf(fid,'& ')
     %fprintf(fid, additions{1,jj});

     fprintf(fid,'\\\\  \n');
    
end

fprintf(fid,' \\cline{1-3}  \n');
for jj=size(Keys,1)-2:size(Keys,1)
    
     fprintf(fid,Keys{jj,1});
     fprintf(fid,' ');
     if values(jj,1)<2000
     fprintf(fid,strcat('& ',num2str(values(jj,1),'%5.3f')));
     else
     fprintf(fid,strcat('& ','$\\infty$'));
     end
     
     %fprintf(fid,strcat('& ',num2str(data_m(jj,1),'%5.3f')));
     fprintf(fid,' ');
     fprintf(fid,strcat('& ',num2str(calibrations(jj,1),'%5.0f')));
     fprintf(fid,' ');
     %fprintf(fid,'& ')
     %fprintf(fid, additions{1,jj});

     fprintf(fid,'\\\\  \n');
    
end
 
%fprintf(fid,'& &&&& &&& &&&  \\\\  \\hline \\hline\n');
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\vspace{0.1in} \\\\  \n');
fprintf(fid,'\\end{center} \n');
%fprintf(fid,'\\justify \n');
fprintf(fid,'\\caption{ \\footnotesize  ');
fprintf(fid,' Baseline parameter values.} \\label{tab:params} ');
 fprintf(fid,'  \n');
fprintf(fid,'\\end{table} \n');
fclose(fid);

cd 'C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Codes\matlab\'


%movefile('MomentsMain.tex','C:\Users\ranos\OneDrive - Umich\Documents\D\Michigan\Res\Female careers in location\Output\tables\');
end  
