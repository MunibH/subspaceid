function rez = subspaceid_psid(meta, obj, params)

%% PREPROCESS
obj.psth = obj.psth(:,:,params.conditions);
rez.time = obj.time;
rez.psth = obj.psth;
rez.trialpsth = obj.trialpsth;
rez.condition = meta.condition(params.conditions);

[obj, meta] = preprocess_maxdiff(meta, obj);




%% VIDEO DATA



end % subspaceid_psid








