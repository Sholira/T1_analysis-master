% This script calculates the T1 curve and value of a set of measurements
% stored in different files in the same folder.
%
% The script loads all the files with extension *.t1 located in the folder
% file_path. Then, it reads all the information in the file and the raw
% data. That info is used to process the curves and calculate the
% spin lattice relaxation time.
%
% Dependencies:
%   T1curve.m
% 
% Changelog:
% 2019.05.09 v0.1 : A new branch start from the scripts used previously to
%                   analysis the data. The intention is to add the curve
%                   calculation, optimization and calculation to this script.
%                   Also, to clean the old set of scripts. This is necesary
%                   to improve the statistical analysis of the data and to
%                   develop the "moving window" analysis. fpm.-
%
% 2019.10.11 v0.2 : It was implemented the reading of binary files. The
%                   script can now read text and binary files and process
%                   the data. Working with binary files is much faster than
%                   using text files, so it is promoted. fpm.-
%
% 2019.10.14 v0.3 : The scripts were moved to a new folder, old files were
%                   removed. The original folder is:
%                   D:\felipe\Microscope\Software\T1_analysis_Matlab
%                   The source code was cleaned.
%                   Source code uploaded to Github. fpm.-
%
% 2021.01.17 v0.4 : The stability of the fiting of a T1 curve was improved
%                   by decreasing the tolerances and increasing the maximum
%                   iterations. Also, for the rolling window calculation,
%                   the result of one fitting is now given as initial 
%                   condition to the fitting of the next window.
%                   modified files: ProcessSegments.m and biexp_model.m

%% Folder to analyse
% Path to the folder where the data files are.
file_path = 'D:\Data\Copper\2020-07-22-copper measurments\measurement 5\combined\H2O2\2020-07-22-FND-H2O2_P01_T1_02\';

%% Parameters
PROCESS_ALL = false; % Activates / deactivates the calculation and the fitting of the T1 curve.
SEGMENT = false; % Activates / deactivates the segmentation analysis of the data.
DYNAMICS = true; % Activates / deactivates the runing window T1 calculation.
MERGE = true; % Activates / deactivates wheter the data should be treated sequentially or not


q_seg = 24; % Quantity of segments to chop the data to.

%% Load the data files. This is done only the first time you run the script.
% If you want to force the loading, delete (clear) the FILES_LOADED variable.
if exist('FILES_LOADED')== 0
    filenames = dir(strcat(file_path,'*.t1')); % get the list of 't1' files in the folder.
    file = strcat(file_path,{filenames.name}); % build the complete file path.
    nf = length(file); % quantity of *.t1 files in the folder.
    % Read all the data files
    for i=1:nf
        T1_experiment(i) = T1curve(cell2mat(file(i)));
    end
    
    % Clear the variables
    clear filenames file nf i
    
    % Avoid to load, again, the files the next time the scrip is run.
    FILES_LOADED = true;
else
    %% ========================== Clear variables =========================
    clear Tcseg Tcseg_o Tcseg T1s;

end
if MERGE == true
    T1_back = T1_experiment;
    timetags = [];
    triggers = [];
    repetitions = 0;
    for i = 1:length(T1_experiment)
        timetags = [timetags; T1_experiment(i).T1_data.timetags];
        triggers = [triggers; T1_experiment(i).T1_data.triggers];
        repetitions = repetitions + T1_experiment(i).T1_data.repetitions;
    end
    T1_new = T1_experiment(1);
    T1_new.T1_data.timetags = timetags;
    T1_new.T1_data.triggers = triggers;
    T1_new.T1_data.repetitions = repetitions;
    T1_new.repetitions = repetitions;
    T1_experiment = T1_new;
end

% Get the numbre of files processed.
nc = length(T1_experiment);

%% ===================== Process all the curves ===========================
if PROCESS_ALL == true
    fprintf('=== Calculating curves ===');
    for i = 1:nc
        T1_experiment(i).process_all;
    end
end
%% ====================== Perform segmentation ============================
if SEGMENT == true
    for i = 1:nc
        Tcseg(i,:) = T1_experiment(i).getSegments(q_seg,true,true);
        txt = sprintf('Segmentation - %s (%d)',T1_experiment(i).name,q_seg);
        title(txt,'Interpreter','none');
    end
end

%% ====================== Calculate dynamics ============================
if DYNAMICS == true
    % Local parameters
    win_size = 6000;
    shift_size = 10;
    
    % Calculate the time of each sample.
    total_time = 20; % minutes
    
    % Quantity of segments
    q_seg = floor(1 + (T1_experiment(1).repetitions-win_size)/shift_size);
    
    T1s = zeros(nc,q_seg);
    error = zeros(nc,q_seg);
    fitparam = zeros(nc,q_seg,5);
    
    for i = 1:nc
        [T1s(i,:), error(i,:), fitparam(i,:,:)] = T1_experiment(i).getDynamics(win_size,shift_size);
        
        % Plot
        t = linspace(0,total_time,length(T1s));
        figure;
        plot(t,T1s(i,:),'.r');
        grid on;
        txt = sprintf('Moving T1 calculation - %s (%d ; %d)',T1_experiment(i).name,win_size,shift_size);
        title(txt,'Interpreter','none');
        xlabel('Time minutes');
        ylabel('T_1 constant [us]');
        lbo = prctile(T1s,25) - (1.5*iqr(T1s));
        ubo = prctile(T1s,75) + (1.5*iqr(T1s));
        T1sno = T1s(T1s > lbo & T1s < ubo);
        tno = t(T1s > lbo & T1s < ubo);
        figure;
        plot(tno,T1sno,'.r')
        grid on;
        txt2 = sprintf ('Moving T1 calculations outliers reloved - %s (%d ; %d)', T1_experiment(i).name, win_size, shift_size);
        title(txt2, 'Interpreter', 'none');
        xlabel('Time minutes');
        ylabel('T_1 constant [us]');
        T = table(t', T1s');
        S = table (tno', T1sno');
        docnam = strcat(file_path,T1_experiment(i).name,'_T1_dynamic.xlsx');
        docnam2 = strcat(file_path,T1_experiment(i).name,'_T1_dynamic_no_outliers.xlsx');
        writetable(T,docnam);
        writetable(S,docnam2);
    end
end
