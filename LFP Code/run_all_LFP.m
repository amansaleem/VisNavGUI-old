function run_all_LFP(run_num, movingwin)

if nargin<2
    movingwin=[3 1];
end

switch run_num
    case 1
        iseries = 1;
        expt_list.animal{1}.name = 'M140501_BALL';
        expt_list.animal{1}.series{iseries}.name = 530;
        expt_list.animal{1}.series{iseries}.expts = [101:108];
        iseries = iseries + 1;
        expt_list.animal{1}.series{iseries}.name = 531;
        expt_list.animal{1}.series{iseries}.expts = [102:109];
        iseries = iseries + 1;
        expt_list.animal{1}.series{iseries}.name = 601;
        expt_list.animal{1}.series{iseries}.expts = [101:109];
        iseries = iseries + 1;
        expt_list.animal{1}.series{iseries}.name = 602;
        expt_list.animal{1}.series{iseries}.expts = [101:109];
        iseries = iseries + 1;
        
    case 2
        iseries = 1;
        expt_list.animal{1}.name = 'M140502_BALL';
        expt_list.animal{1}.series{iseries}.name = 603;
        expt_list.animal{1}.series{iseries}.expts = [106:111];
        iseries = iseries + 1;
        expt_list.animal{1}.series{iseries}.name = 604;
        expt_list.animal{1}.series{iseries}.expts = [106:111];
        iseries = iseries + 1;
        expt_list.animal{1}.series{iseries}.name = 605;
        expt_list.animal{1}.series{iseries}.expts = [106:111 113];
end

for aidx = 1:length(expt_list.animal)
    animal = expt_list.animal{aidx}.name;
    for sidx = 1:length(expt_list.animal{aidx}.series)
        iseries = expt_list.animal{aidx}.series{sidx}.name;
        for eidx = 1:length(expt_list.animal{aidx}.series{sidx}.expts)
            iexp = expt_list.animal{aidx}.series{sidx}.expts(eidx);
            chn_list = 32 + [1 6 10:4:22 28 32];
            es_LFP = VR_LFP_power_safe_allChn(animal,   iseries,   iexp, chn_list, 1, movingwin);
            es_beh = VRLoadMultipleExpts(animal,   iseries,   iexp);
            
            es_LFP.chn_list = chn_list;
            
            saveDir = ['\\zserver.ioo.ucl.ac.uk\Lab\Share\Krumin\Gamma' filesep animal filesep num2str(iseries)];
            fileName = [animal '_' num2str(iseries) '_' num2str(iexp) '_gamma_moreChn_winSize_' num2str(movingwin(1)*1000) '_winShift_' num2str(movingwin(2)*1000)];
            save([saveDir filesep fileName], '-v7.3', 'es_LFP', 'es_beh');
            
            clear es_LFP es_beh;
        end
    end
end
end