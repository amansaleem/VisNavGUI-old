function getAllComodulograms

ianimal = 140511;
iseries = 715;
iexp    = 102;
run_this

ianimal = 140502;
iseries = 605;
iexp    = [107:108]; % g- 113;
run_this

ianimal = 140502;
iseries = 604;
iexp    = [107:108]; % g- 111;
run_this

ianimal = 140502;
iseries = 603;
iexp    = [107:108];
run_this

ianimal = 140501;
iseries = 602;
iexp    = [102:103];
run_this

ianimal = 140501;
iseries = 601;
iexp    = [103:104];
run_this

ianimal = 140501;
iseries = 531;
iexp    = [103:104];
run_this

ianimal = 140501;
iseries = 530;
iexp    = [102:103]; % g- 108;
run_this

ianimal = 130918;
iseries = 1030;
iexp    = [102:103]; % g- 106;
run_this

ianimal = 130920;
iseries = 1025;
iexp    = [102:103]; % g- 106;
run_this

ianimal = 130702;
iseries = 807;
iexp    = [103:104]; % g- 112;
run_this

ianimal = 130703;
iseries = 809;
iexp    = [107:108]; % g- 106;
run_this

    function run_this
        animal = ['M' num2str(ianimal) '_BALL'];
        es = VRLoadMultipleExpts(animal, iseries, iexp,'BEHAV_LFP'); 
        % g- [es] = VR_LFP_power_safe_allChn(animal, iseries, iexp, [10 20 40 50]);
        
        frange = es.freq>1 & es.freq<90;
%g-         corr_CV =  corr(squeeze(es.powB(2,1,:,frange)),squeeze(es.powB(4,1,:,frange)));
%g-         corr_VV =  corr(squeeze(es.powB(4,1,:,frange)),squeeze(es.powB(4,1,:,frange)));
        corr_CV =  corr(es.powA,es.powB);
        corr_VV =  corr(es.powB,es.powB);
        corr_CC =  corr(es.powA,es.powA);
        f = es.freq(frange);
        
        save([animal '_' num2str(iseries) '_corr'], 'corr_CV', 'corr_CC', 'corr_VV', 'f');
    end
end