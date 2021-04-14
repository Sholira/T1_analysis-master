function T1_data = read_bin(filepath)
    %% Show the file to be read.
    fprintf('Reading %s \t|', filepath);

    %% Create a structure array to store the data.
    T1_data.pulse_width = NaN;
    T1_data.read_width = NaN;
    T1_data.read_delay = NaN;
    T1_data.repetitions = NaN;
    T1_data.sampling_mode = NaN;
    T1_data.correction = NaN;
    T1_data.darktimes = [];
    T1_data.PL_measured = [];
    T1_data.filename = '';
    T1_data.timetags = {};
    T1_data.triggers = {};
    T1_data.timestamp = NaN;
    T1_data.elapsed_time = NaN;
    
    %% Get data from file
    % Open the file and get the file name
    [~,T1_data.filename,~] = fileparts(filepath); % get filename.
    fid = fopen(filepath,'r'); % open file to read.
     % Fill one part in the progress bar
    fprintf('-');
    
    % Get the data from the file.
    % Labview saves the files in big-endiand format.
    if fid ~= -1
        % Check that the file is a binary T1 file (0xFE).
        file_type = fread(fid, [1,1],'uint8','ieee-be');
        if file_type == 254 % = 0xFE
            % Fill second part in the progress bar
            fprintf('-');
            
            %% Read the binary file and get the data
            % Build the date. The fraction of seconds are readed but not processed.
            % There is an issue witht the summer-time saving.
            timestp_s = fread(fid, [1,1],'int64','ieee-be'); % seconds after the Epoch
            timestp_fs = fread(fid, [1,1],'uint64','ieee-be'); % fraction of a second
            T1_data.timestamp = datetime(timestp_s,'ConvertFrom','epochtime',...
                                        'Epoch','1904-01-01','TimeZone','local');
   
            T1_data.pulse_width = fread(fid, [1,1],'single','ieee-be');
            T1_data.read_width = fread(fid, [1,1],'single','ieee-be');
            T1_data.read_delay = fread(fid, [1,1],'single','ieee-be');
            T1_data.repetitions = fread(fid, [1,1],'int32','ieee-be');
            T1_data.sampling_mode = fread(fid, [1,1],'uint16','ieee-be');
            T1_data.correction = fread(fid, [1,1],'single','ieee-be');
            
            % Build the elapsed time. The fraction of seconds are readed but not processed.
            timestp_s = fread(fid, [1,1],'int64','ieee-be'); % seconds after the Epoch
            timestp_fs = fread(fid, [1,1],'uint64','ieee-be'); % fraction of a second
            T1_data.elapsed_time = timestp_s + timestp_fs; % in seconds. But check this!!
            
            % Read dark times.
            arraysize = fread(fid, [1,1],'uint32','ieee-be'); % size of the vector containing the darktimes
            T1_data.darktimes = fread(fid, [arraysize,1],'single','ieee-be')';
            
            % Read the T1curve
            arraysize = fread(fid, [1,1],'uint32','ieee-be'); % size of the vector containing the T1 curve
            T1_data.PL_measured = fread(fid, [arraysize,1],'single','ieee-be');
            
            % Read the triggers (2D array)
            colums = fread(fid, [1,1],'uint32','ieee-be'); % columns of the array containing the triggers
            rows = fread(fid, [1,1],'uint32','ieee-be'); % rows of the array containing the triggers
            T1_data.triggers = fread(fid, [rows,colums],'uint32','ieee-be')';
            
            % Read the timetags (2D array)
            colums = fread(fid, [1,1],'uint32','ieee-be'); % columns of the array containing the triggers
            rows = fread(fid, [1,1],'uint32','ieee-be'); % rows of the array containing the triggers
            T1_data.timetags = fread(fid, [rows,colums],'uint32','ieee-be')';
            
                        
            % Fill third part in the progress bar
            fprintf('-');
            
            %% Convert the data types
            T1_data.triggers = num2cell(T1_data.triggers,2);
            T1_data.timetags = num2cell(T1_data.timetags,2);
            
            % Remove zeros from the timetags
            T1_data.timetags = cellfun(@(x) x(x~=0),T1_data.timetags,'UniformOutput',false);
            
            % Fill the last part of the progress bar
            disp('- OK')
        else
            disp('This is not a binary T1 file');
        end
        fclose(fid);
    else
        disp('File reading error');
        return
    end
end