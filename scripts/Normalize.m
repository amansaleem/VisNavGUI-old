function vec_out = Normalize(vec_in)
m = mean(vec_in(:),1);
sd = std(vec_in(:),[],1);
if sd == 0
    sd = 1;
end
vec_out = (vec_in-m)/sd;
end