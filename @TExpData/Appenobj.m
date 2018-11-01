function obj = Appenobj(obj,obj_in)
    %concatenate some of the properties. We shouldn't have to concatenate 
    %the other ones
    obj.Bayes.goodidx = cat(1,obj.Bayes.goodidx,obj_in.Bayes.goodidx);
    obj.Bayes.X = cat(1,obj.Bayes.X,obj_in.Bayes.X); 
    obj.Bayes.X0 = cat(1,obj.Bayes.X0,obj_in.Bayes.X0); 
    for iprobe = 1:size(obj_in.Bayes.Posterior,1)
        for phsidx = 1:size(obj.Bayes.Posterior,2)
            obj.Bayes.Posterior{iprobe,phsidx} = cat(1,obj.Bayes.Posterior{iprobe,phsidx},obj_in.Bayes.Posterior{iprobe,phsidx});
            obj.Bayes.PosError{iprobe,phsidx} = cat(1,obj.Bayes.PosError{iprobe,phsidx},obj_in.Bayes.PosError{iprobe,phsidx});
            obj.Bayes.DistError{iprobe,phsidx} = cat(1,obj.Bayes.DistError{iprobe,phsidx},obj_in.Bayes.DistError{iprobe,phsidx});
            obj.Bayes.prediction{iprobe,phsidx} = cat(1,obj.Bayes.prediction{iprobe,phsidx},obj_in.Bayes.prediction{iprobe,phsidx});
        end
        obj.Bayes.Posterior0{iprobe} = cat(1,obj.Bayes.Posterior0{iprobe},obj_in.Bayes.Posterior0{iprobe});
        obj.Bayes.PosError0{iprobe} = cat(1,obj.Bayes.PosError0{iprobe},obj_in.Bayes.PosError0{iprobe});
        obj.Bayes.DistError0{iprobe} = cat(1,obj.Bayes.DistError0{iprobe},obj_in.Bayes.DistError0{iprobe});
        obj.Bayes.LFPphase{iprobe} = cat(1,obj.Bayes.LFPphase{iprobe},obj_in.Bayes.LFPphase{iprobe});
    end
    obj.Bayes.D = cat(1,obj.Bayes.D,obj_in.Bayes.D);    
    obj.Bayes.time = cat(1,obj.Bayes.time,obj.Bayes.time(end) + obj_in.Bayes.time);    
    obj.Bayes.Kroom = cat(1,obj.Bayes.Kroom,obj_in.Bayes.Kroom);
    obj.Bayes.Kdist = cat(1,obj.Bayes.Kdist,obj_in.Bayes.Kdist);
    
    obj.Bayes.kalVisstd = cat(1,obj.Bayes.kalVisstd,obj_in.Bayes.kalVisstd);
    
    obj.data.es = combineTwoVRseries(obj.data.es, obj_in.data.es);
end