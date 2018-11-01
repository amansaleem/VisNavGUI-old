function time = convertdate2time(date)
time = date(:,end-2)*3600 + date(:,end-1)*60 + date(:,end);
end