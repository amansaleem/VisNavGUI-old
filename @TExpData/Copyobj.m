function obj = Copyobj(obj,obj_in)
    prop = properties(obj_in);
    for p = 1:numel(prop)
        obj.(prop{p})=obj_in.(prop{p});    
    end
end