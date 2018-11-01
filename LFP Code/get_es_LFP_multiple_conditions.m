parfor n = 12:13
    switch n
        case 1
            VRLoadMultipleExpts('M130920_BALL',1025,[106 108],'BEHAV_LFP','gray');
        case 2
            es_OL = VRLoadMultipleExpts('M130920_BALL',1025,[201:203],'BEHAV_LFP','OL');
        case 3
            es_CL1 = VRLoadMultipleExpts('M130920_BALL',1025,[101:105],'BEHAV_LFP','CL1');
        case 4
            es_CL2 = VRLoadMultipleExpts('M130920_BALL',1025,[107],'BEHAV_LFP','CL2');
        case 5
            es2_gray = VRLoadMultipleExpts('M130918_BALL',1030,[106 109],'BEHAV_LFP','gray');
        case 6
            es2_OL = VRLoadMultipleExpts('M130918_BALL',1030,[201:202],'BEHAV_LFP','OL');
        case 7
            es2_CL1 = VRLoadMultipleExpts('M130918_BALL',1030,[101:105],'BEHAV_LFP','CL1');
        case 7
            es2_CL2 = VRLoadMultipleExpts('M130918_BALL',1030,[107 108 110],'BEHAV_LFP','CL2');
        case 9
            VRLoadMultipleExpts('M130703_BALL',809,101:104,'BEHAV_LFP','CL1');
        case 10
            VRLoadMultipleExpts('M130703_BALL',809,106:109,'BEHAV_LFP','CL2');
        case 11
            VRLoadMultipleExpts('M130703_BALL',809,110:113,'BEHAV_LFP','CL3');
        case 12
            VRLoadMultipleExpts('M130702_BALL',807,102:107,'BEHAV_LFP','CL1');
        case 13
            VRLoadMultipleExpts('M130702_BALL',807,109:112,'BEHAV_LFP','CL2');
    end
end