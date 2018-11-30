%% Tomaso Muzzu - UCL _ 28/02/2018

% function to save the raw data in a dat file


function save2dat(EphysData, rec_nr, FileName)

if rec_nr == 1
    % Save data    
    fid = fopen(FileName, 'w');
    fwrite(fid, EphysData,'int16');
    fclose(fid);
else     
    % Save data
    fid = fopen(FileName, 'a');
    fwrite(fid, EphysData,'int16'); 
    fclose(fid); 
end

end