classdef TExpData < handle
%    
    properties
        animal
        series
        iseries
        exptList
        data
        CellInfo
    end
    
    methods
        function obj = TExpData()
            obj.animal = [];
            obj.series = [];
            obj.iseries = [];
            obj.exptList = [];
            obj.data = [];
            obj.CellInfo = [];
        end
        
        obj = SelectAnimal(obj, DIRS, strcell);            
        obj = SelectExpt(obj, Fall);
        obj = SelectSeries(obj, selmode);
        iseries_list = getallSeries(obj);
        obj = defineCellInfo(obj, spikeIDs, chanIDs, probeIDs);
        obj = Copyobj(obj,obj_in);
        obj = Appenobj(obj,obj_in);
    end
end
