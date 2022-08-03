
function params=untransform(params)

% untransform
params('LA0','value')={1/(exp(-params{'LA0',:}) +1)};
params('plocal_','value')={1/(exp(-params{'plocal_',:}) +1)};
params('piw_','value')={1/(exp(-params{'piw_',:}) +1)};
%params('Jseg','value')={1/(exp(-params{'Jseg',:}) +1)};

params('deltaw_','value')={1/(exp(-params{'deltaw_',:}/10) +1)};

params('line','value')={(1+ 1/(exp(-(params{'line',:}-1)) +1))  };


params('pish_','value')={exp(params{'pish_',:})};
params('PHID_','value')={exp(params{'PHID_',:})};

params('eta_','value')={exp(params{'eta_',:})};
params('ce_','value')={exp(params{'ce_',:})};
params('ces_','value')={exp(params{'ces_',:})}; % irrelevant
params('le_','value')={exp(params{'le_',:})};

params('sel','value')={exp(params{'sel',:})};
params('Jcenter','value')={exp(params{'Jcenter',:})};
params('jobdif_','value')={exp(params{'jobdif_',:})};
params('muw','value')={exp(params{'muw',:})};
params('sw','value')={exp(params{'sw',:})};
params('Xi','value')={exp(params{'Xi',:})};
params('crra_','value')={exp(params{'crra_',:})};
params('crrat_','value')={exp(params{'crrat_',:})};
params('pid_','value')={exp(params{'pid_',:})};
params('wc','value')={exp(params{'wc',:})};
params('crrah_','value')={exp(params{'crrah_',:})};
params('crrax_','value')={exp(params{'crrax_',:})};
params('piel_','value')={exp(params{'piel_',:})};
params('pi_','value')={exp(params{'pi_',:})};
params('NLY','value')={exp(params{'NLY',:})};


params('wgap_raw','value')={-exp(params{'wgap_raw',:})};


%params('pi_','value')={exp(params{'pi_',:})};


params('d','value')={exp(params{'d',:})*10};
%params('gamd_','value')={exp(params{'gamd_',:})*100};
end