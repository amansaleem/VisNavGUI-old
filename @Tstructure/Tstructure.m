classdef Tstructure < dynamicprops
%
    properties
        name
    end
    
    methods
        function obj = Tstructure(str)
            obj.name = str;
        end
        
        %inherited method from dynamicprops:
        %obj = addprop(obj,'PropName')
    end
end
