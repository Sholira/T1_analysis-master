function T1_data = read_T1_file(filepath)
    %% Show the file to be read.
    fprintf('Reading %s \t|', filepath);
    
    %% Script's parameters
    hdr_size = 15; % Quantity of lines of the header. Include the vectors darktimes, PL_measured, PL_fit.

    %% Create a structure array to store the data.
    T1_data.pulse_width = NaN;
    T1_data.read_width = NaN;
    T1_data.read_delay = NaN;
    T1_data.repetitions = NaN;
    T1_data.sampling_mode = NaN;
    T1_data.correction = NaN;
    T1_data.darktimes = [];
    T1_data.PL_measured = [];
    T1_data.PL_fit = [];
    T1_data.filename = '';
    T1_data.timetags = {};
    T1_data.triggers = {};
    
    %% Get data from file
    % Open the file and get the file name.
    [~,T1_data.filename,~] = fileparts(filepath); % get filename.
    fid = fopen(filepath,'rt'); % open file to read, in text mode.
    
    % Read the complete file as an array of chars.
    strData = fscanf(fid,'%c');
    fclose(fid);
    nCharPerLine = diff([0 find(strData == newline) numel(strData)]); % calculate the length of each line.
    cellData = strtrim(mat2cell(strData,1,nCharPerLine)); % get each line as a cell in an array.

    % Delete the unused variables
    clear strData nCharPerLine
    
    %% Fill the one part of the progress bar
    fprintf('-');
    
    %% Get the value of the parameters
    % Parse line by line of the header to find the values and convert them to
    % numbers. Assign the value to the corresponding parameter.
    T1_data.pulse_width = get_value(cell2mat(cellData(1)));
    T1_data.read_width = get_value(cell2mat(cellData(2)));
    T1_data.read_delay = get_value(cell2mat(cellData(3)));
    T1_data.repetitions = get_value(cell2mat(cellData(4)));
    T1_data.sampling_mode = get_value(cell2mat(cellData(5)));
    T1_data.correction = get_value(cell2mat(cellData(6)));
    % The remaining parameters in the header are not used anymore.

    %% Get the saved ralaxation curve
    % Extract the darktimes, the previously calculated curve and fit.
    T1_data.darktimes = str2double(strsplit(cell2mat(cellData(13))));
    T1_data.PL_measured = str2double(strsplit(cell2mat(cellData(14))));
    T1_data.PL_fit = str2double(strsplit(cell2mat(cellData(15))));

    %% Fill the one part of the progress bar
    fprintf('-');
    
    %% Get the raw data
    % The next 'repetitions' lines store the values of the triggers.
    start_trig = hdr_size + 1;
    end_trig = T1_data.repetitions + hdr_size;
    T1_data.triggers = cellfun(@cellstr2vec,cellData(start_trig:end_trig),'UniformOutput',false)';

    % Fill the one part of the progressing bar
    fprintf('-');
    
    % The timetags are the last block of data in the file.
    end_tt = end_trig + T1_data.repetitions;
    T1_data.timetags = cellfun(@cellstr2vec,cellData(end_trig+1:end_tt),'UniformOutput',false)';
    
    %% Fill the last part of the progress bar
    disp('- OK')
end