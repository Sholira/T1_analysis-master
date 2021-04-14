function T1_data = read_T1_file(filepath)
    %% Check if the file is a text or binary file
    % Open the file
    fid = fopen(filepath,'r'); % open file to read.
    
    % Read the first byte
    file_type = fread(fid, [1,1],'uint8','ieee-be');
    
    % Close the file
    fclose(fid);
    
    % If the byte value is 254, the file is binary. If the value is 80,
    % then the file contains text.
    switch file_type
        case 254
            T1_data = read_bin(filepath);
            % Fill the missing fields with default data
            T1_data.PL_fit = [];
       
        case 80
            T1_data = read_txt(filepath);
        
        otherwise
            T1_data = {};
            disp('Wrong T1 file');
    end
end