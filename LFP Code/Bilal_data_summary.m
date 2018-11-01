% %% Across Recordings significance
% % EPSCs and coupled LFP
[fr, residual_power_el, peak_pos_el, peak_val_el, p_el]= narrowgamma_eval(f, sq(nanmean(all_el_spont,3)), [55 70], [20 90]);
[fr, residual_power_ep, peak_pos_ep, peak_val_ep, p_ep]= narrowgamma_eval(f, sq(nanmean(all_ep_spont,3)), [55 70], [20 90]);
% % IPSCs and coupled LFP
[fr, residual_power_il, peak_pos_il, peak_val_il, p_il]= narrowgamma_eval(f, sq(nanmean(all_il_spont,3)), [55 70], [20 90]);
[fr, residual_power_ip, peak_pos_ip, peak_val_ip, p_ip]= narrowgamma_eval(f, sq(nanmean(all_ip_spont,3)), [55 70], [20 90]);
%% Significance for each recording
% And difference in peak positions
for iRec = 1:size(all_ip_spont,3)
    [~, residual_el(iRec,:,:), peak_freq_el(iRec,:), peak_height_el(iRec,:), ps_el(iRec)]= narrowgamma_eval(f, squeeze(all_el_spont(:,:,iRec)), [55 75], [20 90], [19 20]);
    [~, residual_ep(iRec,:,:), peak_freq_ep(iRec,:), peak_height_ep(iRec,:), ps_ep(iRec)]= narrowgamma_eval(f, squeeze(all_ep_spont(:,:,iRec)), [55 75], [20 90], [19 20]);
    %
    [~, residual_il(iRec,:,:), peak_freq_il(iRec,:), peak_height_il(iRec,:), ps_il(iRec)]= narrowgamma_eval(f, squeeze(all_il_spont(:,:,iRec)), [55 75], [20 90], [19 20]);
    [~, residual_ip(iRec,:,:), peak_freq_ip(iRec,:), peak_height_ip(iRec,:), ps_ip(iRec)]= narrowgamma_eval(f, squeeze(all_ip_spont(:,:,iRec)), [55 75], [20 90], [19 20]);
end
    
%% For Figure 1D
es_11 = VR_LFP_power_new('M140212_BALL', 222, 11, 10, 0, [3 1], []);
es_12 = VR_LFP_power_new('M140212_BALL', 222, 12, 10, 0, [3 1], []);
es_13 = VR_LFP_power_new('M140212_BALL', 222, 13, 10, 0, [3 1], []);
es_14 = VR_LFP_power_new('M140212_BALL', 222, 14, 10, 0, [3 1], []);
es_16 = VR_LFP_power_new('M140212_BALL', 222, 16, 10, 0, [3 1], []);
es_103 = VR_LFP_power_new('M140212_BALL', 222, 103, 10, 0, [3 1], []);
es_102 = VR_LFP_power_new('M140212_BALL', 222, 102, 10, 0, [3 1], []);

[fr, residual_11, peak_pos_11, peak_val_11, p_11]= narrowgamma_eval(es_11.freq, es_11.powA, [55 70], [20 90]);
[fr, residual_12, peak_pos_12, peak_val_12, p_12]= narrowgamma_eval(es_12.freq, es_12.powA, [55 70], [20 90]);
[fr, residual_13, peak_pos_13, peak_val_13, p_13]= narrowgamma_eval(es_13.freq, es_13.powA, [55 70], [20 90]);
[fr, residual_14, peak_pos_14, peak_val_14, p_14]= narrowgamma_eval(es_14.freq, es_14.powA, [55 70], [20 90]);
[fr, residual_16, peak_pos_16, peak_val_16, p_16]= narrowgamma_eval(es_16.freq, es_16.powA, [55 70], [20 90]);
[fr, residual_102, peak_pos_102, peak_val_102, p_102]= narrowgamma_eval(es_102.freq, es_102.powA, [55 70], [20 90]);
[fr, residual_103, peak_pos_103, peak_val_103, p_103]= narrowgamma_eval(es_103.freq, es_103.powA, [55 70], [20 90]);