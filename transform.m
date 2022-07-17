
function params=transform(params)

% untransform
flogit=@(x) -log(1/x -1);
params('LA0','value')={flogit(params{'LA0',:})};
params('plocal_','value')={flogit(params{'plocal_',:})};
params('piw_','value')={flogit(params{'piw_',:})};
%params('Jseg','value')={flogit(params{'Jseg',:})};

params('deltaw_','value')={flogit(params{'deltaw_',:})*1000};


params('line','value')={ (1+ flogit((params{'line',:}-1))) };
params('pish_','value')={log(params{'pish_',:})};

params('PHID_','value')={log(params{'PHID_',:})};
params('eta_','value')={log(params{'eta_',:})};
params('ce_','value')={log(params{'ce_',:})};
params('ces_','value')={log(params{'ces_',:})}; 

params('sel','value')={log(params{'sel',:})};
params('Jcenter','value')={log(params{'Jcenter',:})};
params('jobdif_','value')={log(params{'jobdif_',:})};
params('muw','value')={log(params{'muw',:})};
params('sw','value')={log(params{'sw',:})};
params('Xi','value')={log(params{'Xi',:})};
params('crra_','value')={log(params{'crra_',:})};
params('crrat_','value')={log(params{'crrat_',:})};
params('pid_','value')={log(params{'pid_',:})};
params('wc','value')={log(params{'wc',:})};
params('crrah_','value')={log(params{'crrah_',:})};
params('crrax_','value')={log(params{'crrax_',:})};
params('piel_','value')={log(params{'piel_',:})};
params('pi_','value')={log(params{'pi_',:})};
params('NLY','value')={log(params{'NLY',:})};


params('wgap_raw','value')={ log((-params{'wgap_raw',:}))};

params('d','value')={log(params{'d',:}/10)};
%params('gamd_','value')={log(params{'gamd_',:}/100)};
end