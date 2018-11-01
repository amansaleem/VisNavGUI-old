function obj = CalculateStimTuning(obj, exptlist, shank_list, suffix)
[obj.data.VStuning, obj.data.VSstim, ~, ~] = getStimTuning(obj.animal, obj.iseries, exptlist, shank_list, suffix);
end