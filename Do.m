%{
Mean distance to 'opportunities'/'jobs' given a selected distribution of
jobs and a location

%}

function Do=Do(i,dist,D)

Do=sum(D(:,i).*dist) ; %number

end