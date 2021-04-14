clear all;clc;
rng('shuffle')


%% Folder to analyse
% Path to the folder where the data files are.
%file_path = 'D:\Data\Petri coat\2021-04-02_petricoat_FND-FN\combined\';
file_path = 'D:\Data\Data analysis\Data from gadolinium\bootstrap repeat 1 2500\';
%file_path = 'D:\Data\Polymer\2021-03-09_polymer\combined\';
%file_path = 'D:\Data\Petri coat\2021-03-09_petricoat_FND-HA\combined\';

%% Parameters

MODEL       = 'all'; % biexp (biexponential model) | siexp (single exponential) | srexp (tretched exponential | all (fit with all models) | false (don't do fitting) || which model will be used during the analysis
FIT         = false; % TRUE or FALSE; turns on or off whether all the data has to be fitted with the specified models
NORM        = true; % TRUE or FALSE; turns on or off for the normalization; stretched exponential will always be normalised
MERGE       = false; % TRUE or FALSE; merge the T1 files to be considered as a single file
SEGM        = false; % TRUE or FALSE; segment the data into a preset number of segments
RoWi        = false; % TRUE or FALSE; apply the rolling window with winsize window size and shift as the shift
BOOTS       = true; % TRUE or FALSE; do a bootstrap of the measurement
DATABASE    = false; % TRUE or FALSE; turns on or off whether a full database will be kept of the bootstrap (the database will be generated but when done, it will be removed)
COMBI       = false; % TRUE or FALSE; do a rolling window combined with bootstrap

%% Variables

bootsam     = 2500; % The number of bootstrap samples
segm        = 2; % The number of segments
window      = 5000; % The number of repetitions used for one T1 in the rolling window
shift       = 100; % The number of repetitions used for the shift of the rolling window
time        = 20; % the total time of the measurement (also of the merge)

%% loading the files

if exist('FILES_LOADED')==0
    filenames = dir(strcat(file_path,'*.t1')); % get the list of the t1 files in the folder
    file = strcat (file_path, {filenames.name}); % build the complete file path
    nf = length (file); %quantity of the T1 files in the folder
    for i = 1:nf
        T1_experiment(i) = T1curve(cell2mat(file(i)));
    end
    clear file nf i % clear the variables
    FILES_LOADED = true; % avoid to load the files again if the script is ran again
else
    clear Tcseg Tcseg_o T1s;
end

%% Fit the models

if FIT == true
    model_fit(T1_experiment, file_path, MODEL, NORM);
end

%% Merging of the data

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

%% Segmentation of the data

if SEGM == true
    T1_segmentation(T1_experiment, segm, MODEL, file_path,NORM);
end

%% Apply the rolling window analysis to the data

if RoWi == true
    T1_rollingwindow(T1_experiment, file_path, MODEL, window, shift, time, NORM, COMBI, bootsam);
end

%% Do a bootstrap of the data

if BOOTS==true
    for i = 1:length(T1_experiment)
        name = T1_experiment(i).name;
        timetags = T1_experiment(i).T1_data.timetags;
        triggers = T1_experiment(i).T1_data.triggers;
        darktimes = T1_experiment(i).darktimes;
        rwid = T1_experiment(i).read_width;
        corr = T1_experiment(i).correction;
        T1_bootstrap(name, timetags, triggers, darktimes, rwid, corr, MODEL, NORM, bootsam, file_path);
    end
end
