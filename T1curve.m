% This class contains the data from a T1 file and defines the methods that
% can operate over it.
% Method:
%   T1curve(filepath) : Constructor, creates the object by reading a T1
%                       file.
%   calc_T1curve()    : Calculates the T1 relaxation curve with the
%                       information stored in the 'T1curve' object.
%   fit_model(model,raw) : Find the curve that best fits the model defined
%                          by the input parameter 'model' (only 'biexponential'
%                          is working now).
%                          model: 'biexponential' (default)
%                                 'exponential'   (not implemented yet)
%                                 'singleNV'      (not implemented yet)
%                          raw: {true, false} if false it uses a normalised
%                          version of the relaxation curve. If true, uses
%                          the curve from the raw data.
%   plot_curve(type): Plots the relaxation curve, the fitting curve and
%                     shows the T1 value in the figure.
%                     type: 'normalised' shows the normalised curve.
%                           'raw' shows the curve from the raw data (default).
%   getPLsignal(fromto) : Integrates all the runs into one sequence of
%                         pulses. The output 'sum_pl(i)' is the addition of
%                         all the photons detected at time 'i' (referenced
%                         to the begining of the current run).
%                         fromto: a vector [a b] defining the initial (a) and
%                                 last (b) run to process. Default, all the
%                                 runs.
%   getPulses(obj,fromto) : Chops the train of pulses produced by
%                           getPLsignal into isolated pulses.
%                           fromto: a vector [a b] defining the initial (a) and
%                                   last (b) run to process. Default, all the
%                                   runs.
%   process_all(obj,plt) : Run, sequentially the methods calc_T1curve and
%                          fit_model. if 'plt' (plot) is 'true' then it
%                          also runs plot_curve.
%   getSegments(segments,norm,plt) : Calculate the T1 constants spliting the 
%                                    dataset to 'segments' subsets.
%
%   getDynamics(win_size,shift_size) : calculates the T1 values over a subset 
%                                      of T1 runs, using the runing window
%                                      process.
%
%   evalQArun(plt) : Evaluate the quality of each run
%
%   reset_filter() : Reset the filtering of the data
%
%   reset_curve() : Reset the curve and the calculations.
%
%   getFFT(signal) : Plot the FFT of the integrated PL.
%
%   getname() : Print the name of all the files in the dataset.
%
%   getdarktimes) : Return the vector containing the darktimes used in the experiment.
%
%   setWhitelist(list) : Sets the runs to be used to calculate the T1 curve.
%
% Dependencies:
%   read_T1_file.m
%   doNORM.m
%   getT1curve.m
%   biexp_model.m
%   processSegmets.m
%
classdef T1curve < handle
    properties
        correction = NaN; % Value of the correction used to syncronise the triggers [us].
        read_width = NaN; % Size of the detection window [us].        
        model = "biexponential"; % Model used to fit the data and calculatethe T1 value.
        name = '' % String. Name to show in legends and tables
        group = '' % String. Name of the data set (group) which the curve belongs.
        exp_num = 1; % Double, number of the measurment in a set of experiments.
        whitelist = []; % Set of runs to use in the analysis.
        QA_score = NaN;        
    end
    
    properties
        darktimes = []; % Corresponding dark times of the T1_curve [us].
        T1_curve = []; % Integrated PL over the detection window of each reading pulse. Curve not normalised
        T1_curve_n = []; % Normalised T1 curve.
        T1_curve_f = []; % Fitted curve acording to the model used.
        T1 = NaN; % Value of the spin - lattice relaxation constant.
        PL_T1 = NaN; % Photoluminescence at T1.
        fit_param = []; % Array with all the calculated parameters used to fit the model.
        fit_error = NaN; % Norm of the residues of the fitting.
        repetitions = NaN; % Amount of repetitions of a single run (one scanning of all the darktimes).
        pulse_width = NaN; % Size of the laser pulse in [us].
        t_int = []; % Vector with the time line of the integrated measurements.
        triggers_int = []; % Calculated triggers (mean value) for the integrated PL signal.
        PL_int = []; % Integrated PL signal, built from all the runs recorded in the experiment.
        pulses_int = []; % Isolated pulses from the integrated PL, each row has a pulse of pulse_width length.
        SNR = NaN; % Signal to Noise Ratio. (From the integrated PL signal).

        T1_data
    end
    
%     properties (Hidden, GetAccess = private, SetAccess=private)
%         T1_data
%     end
    
    
    methods
        %% Constructor
        function obj=T1curve(filepath)
            % Read the file and get the data as a cell array.
            obj.T1_data =  read_T1_file(filepath);
            
            % Fill the properties with the data.
            obj.read_width = obj.T1_data.read_width;
            obj.correction = obj.T1_data.correction + obj.T1_data.read_delay;
            obj.darktimes = obj.T1_data.darktimes;
            obj.repetitions = obj.T1_data.repetitions;
            obj.pulse_width = obj.T1_data.pulse_width;
            obj.T1_curve = obj.T1_data.PL_measured;
            obj.T1_curve_n = obj.T1_data.PL_measured;
            obj.T1_curve_f = obj.T1_data.PL_fit;
            obj.name = obj.T1_data.filename;
            obj.whitelist = 1:obj.T1_data.repetitions;
            tmp = strsplit(obj.T1_data.filename,'_');
            if length(tmp)>1
                obj.exp_num = str2double(tmp(end));
                obj.group = obj.T1_data.filename(1:end-3);
            else
                obj.exp_num = 1;
                obj.group = obj.T1_data.filename;
            end
            
            % Release memory
            clear T1_data tmp
        end
        
        
        %% Calculate the relaxation curve 
        function calc_T1curve(obj)
            % Calculate the raw curve
            data_f = obj.T1_data.timetags(obj.whitelist);
            trigg_f = obj.T1_data.triggers(obj.whitelist);
            obj.T1_curve = getT1curve(trigg_f,data_f,obj.correction,obj.read_width);
            
            % Normalyse the curve
            obj.T1_curve_n = doNORM(obj.T1_curve,'tail');
        end
        
        
        %% Apply the model to calculate the T1 value
        function fit_model(obj,model,raw)
            curve = obj.T1_curve_n;
            if exist('raw','var') ~= 0
                if raw
                    curve = obj.T1_curve;
                end
            end
            
            if exist('model','var') == 0
                model = 'biexponential';
            end

            obj.model = model;
            
            switch model
                case 'biexponential'
                    [obj.T1_curve_f, obj.fit_param, obj.PL_T1, obj.fit_error] = biexp_model(obj.darktimes,curve);
                    obj.T1 = max(obj.fit_param([3 5])); % T1 is the longest time constant.
                case 'exponential'
                    disp('Not implemented yet.')
%                     res_tab_1exp = FIToneEXP_model(obj.darktimes,obj.T1_curve,linename);
                case 'singleNV'
                    disp('Not implemented yet.')
                otherwise
                    disp('Select between biexponential, exponential or singleNV')
            end
        end
        
        
        %% Plot the T1 and fitting curves.
        function plot_curve(obj,type)
            figA = figure;

            % Get the current colour palette (line colors)
            co = get(groot,'DefaultAxesColorOrder');

            % Set the same color to the data points and the fitting curve.
            co = reshape([co(:) co(:)]',2*size(co,1), []); % Assign the same color to two consecutive lines
            set(figA,'defaultAxesColorOrder',co); % Set the new line color array to the axes.

            % Plot the lines
            if exist('type','var') == 0
                type = 'normalised';
            end
            switch type
                case 'raw'
                    semilogx(obj.darktimes,obj.T1_curve,'.');
                    legend(obj.name);
                otherwise
                    semilogx(obj.darktimes,obj.T1_curve_n,'.',obj.darktimes,obj.T1_curve_f,'-',obj.T1,obj.PL_T1,'pr');
                    val = sprintf('T1: %.2f [us]', obj.T1);
                    legend({obj.name, 'Fitted curve', val});
            end
            l=legend('show');
            set(l, 'Interpreter', 'none');
            legend boxoff;
            grid on;
            title('Spin-lattice relaxation curve');
            xlabel('Dark time [us]');
            ylabel('Intensity [a.u.]');
            
        end
        
        
        %% Build the pulsed PL signal from a set of measures
        function [t, triggers, sum_pl] = getPLsignal(obj,fromto)
            % Check the optional parameter
            if exist('fromto','var') == 0
                fromto = [1 obj.repetitions];
            end
            
            % Build the time variable
            t = 1:max(cellfun(@max,obj.T1_data.timetags(fromto(1):fromto(2),:))); % [counter cycles (10ns)]
            
            % Count the amount of photons detected at the same time after started the experiment.
            sum_pl = zeros(1,length(t));
            for i=fromto(1):fromto(2)
                sum_pl(obj.T1_data.timetags{i}) = sum_pl(obj.T1_data.timetags{i}) + 1;
            end
            
            % Average the trigger position
            triggers = [1 mean(cell2mat(obj.T1_data.triggers))] + obj.correction*100; % [index = 10 ns]
            triggers = round(triggers);
            
            % Store the result in the object
            obj.t_int = t;
            obj.triggers_int = triggers;
            obj.PL_int = sum_pl;
        end
        
        
        %% Separate the pulses from getPLsignal.
        function pulses = getPulses(obj,fromto)
            % Check the optional parameter
            if exist('fromto','var') == 0
                fromto = [1 obj.repetitions];
            end
            
            % Amount of pulses in one run
            q_pulses = length(obj.darktimes);
            
            % Convert the pulse width [us] to index
            p_width = obj.pulse_width * 100; % samples
            
            % Get the train of pulses
            [~, triggers, sum_pl] = getPLsignal(obj,fromto);
            
            pulses = zeros(q_pulses,p_width);
            for i = 1:q_pulses
                pulses(i,:) = sum_pl(triggers(i):(triggers(i) + p_width -1));
            end
            
            % Store the result in the object
            obj.pulses_int = pulses;
        end
        
        
        %% Calculate the curve, fit the model and plot them in one shoot.
        function process_all(obj,plt)
            % Check the optional parameter
            if exist('plt','var') == 0
                plt = false;
            end            

            q_files = length(obj);
            for i=1:q_files
                obj(i).calc_T1curve
                obj(i).fit_model;
                if plt
                    obj(i).plot_curve;
                end 
            end
        end
        
        %% Calculate the T1 constants spliting the dataset into 'segments' subsets.
        function T1s = getSegments(obj,segments,norm,plt)
            % Check the optional parameter
            if exist('norm','var') == 0
                norm = false;
            end
            if exist('plt','var') == 0
                plt = false;
            end
            
            % Check if the amount of repetitions can be divided by the
            % segment size
            q_run_seg = floor(obj.T1_data.repetitions/segments);               
            remanent = obj.T1_data.repetitions - q_run_seg*segments;
            if remanent == 0
                % Set the amount of runs in each segment
                q_run_seg = q_run_seg * ones(1,segments);
            else
                % Set the amount of runs in each segment
                q_run_seg = [q_run_seg * ones(1,segments-1), q_run_seg + remanent];
            end
            
            % Find the position at which each segment starts and stops.
            seg_start = zeros(1,segments);
            seg_stop = zeros(1,segments);
            seg_stop_ = 0;
            for i=1:segments
                seg_start(i) = seg_stop_ + 1;
                seg_stop_ = sum(q_run_seg(1:i));
                seg_stop(i) = seg_stop_;
            end
            fromto(:,1) = seg_start;
            fromto(:,2) = seg_stop;
            
            % Segment the data and calculate the curve.
            curve = zeros(segments,length(obj.darktimes));
            for i=1:segments
                seg_timetags = obj.T1_data.timetags(seg_start(i):seg_stop(i));
                seg_triggers = obj.T1_data.triggers(seg_start(i):seg_stop(i));
                tmp = getT1curve(seg_triggers,seg_timetags,obj.correction,obj.read_width);
                curve(i,:) = tmp;
            end
            
            % Calculate the T1 value of each curve
            T1curve_f = zeros(size(curve));
            fitparam = zeros(segments, 5);
            T1s = zeros(segments,1);
            PLT1 = zeros(segments,1);
            fiterror = zeros(segments,1);
            for i=1:segments 
                [T1curve_f(i,:), fitparam(i,:), PLT1(i), fiterror(i)] = biexp_model(obj.darktimes,curve(i,:));
                T1s(i) = max(fitparam(i,[3 5])); % T1 is the longest time constant.
            end
            
            %% Normalyse before ploting if norm==true.
            if norm
                norm_vals = T1curve_f(:,end);
                curve_n = curve./norm_vals;
                curve_f_n = T1curve_f./norm_vals;
                curve = curve_n;
                T1curve_f = curve_f_n;
                PLT1 = PLT1./norm_vals;
            end
            
            %% Show the time constants
            Tc_tab = table((1:segments)',fromto,T1s,'VariableNames',{'Segment','From_To','Tc'});
            disp(Tc_tab);
%             disp(fitparam);
%             disp(fiterror);
            
            %% Plot all the curves
            if plt
                fig_seg = figure;

                % Get the current colour palette (line colors)
                co = get(groot,'DefaultAxesColorOrder');

                % Flip the yellow color with the green
                co([3 5],:) = co([5 3],:);

                % Set the same color to the data points and the fitting curve.
                co = reshape([co(:) co(:)]',2*size(co,1), []); % Assign the same color to two consecutive lines.
                set(fig_seg,'defaultAxesColorOrder',co); % Set the new line color array to the axes.

                for i=1:segments
                    semilogx(obj.darktimes,curve(i,:),'.',obj.darktimes,T1curve_f(i,:),'-',T1s(i),PLT1(i),'pr');
                    % Set the filename to show in the legend.
                    text = sprintf('Segment %d',i);
                    fig_seg.Children.Children(3).DisplayName = text;
                    % Show only the data line in the legend.
                    fig_seg.Children.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
                    fig_seg.Children.Children(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Increase the size of the point.
                    fig_seg.Children.Children(3).MarkerSize = 12;
                    hold on
                end
                grid on;
                title('Spin-lattice relaxation curve');
                xlabel('Dark time [us]');
                ylabel('Intensity [a.u.]');            
                legend boxoff;

                % Plot the constants in a separate figure
                pos_x = mean(fromto,2);
                fig_t1s = figure;
                bar(pos_x,T1s,1);
                xlim([0 obj.T1_data.repetitions])
%                 ylim([0 120]);
                grid on;
                title('Segments T_1 values');
                xlabel('Run #');
                ylabel('T_1 constant [us]');
            end
        end
        
        
        %% getDynamics : calculates the T1 values over a subset of T1 runs.
        %  The inputs are:
        %   win_size: quantity of runs from which to calculate the constant.
        %   shift_size: The amount of runs to move forward to calculate a
        %   T1 constant (step size).
        function [T1s, error, fitparam] = getDynamics(obj,win_size,shift_size)
            % Find the position at which each segment starts and stops.
            fromto(:,1) = 1:shift_size:(obj.T1_data.repetitions - win_size + 1);
            fromto(:,2) = win_size:shift_size:obj.T1_data.repetitions;
            
            [curves, T1s, PLT1, error, fitparam] = processSegmets(obj.T1_data.timetags,obj.T1_data.triggers,obj.darktimes,fromto,obj.correction,obj.read_width);
            
%             % Sort the paramenters, so the longer time constant is in the
%             % fifth position.
%             long = fitparam(:,3) > fitparam(:,5);
%             tmp = fitparam(long,[2 3]);
%             fitparam(long,[2 3]) = fitparam(long,[4 5]);
%             fitparam(long,[4 5]) = tmp;
           
%             curves_n = curves./min(curves,[],2);
%             plot(obj.darktimes,curves_n);
        end
        
        
        %% Evaluate the quality of each run
        function score = evalQArun(obj,plt)
            if exist('plt','var') == 0
                plt = false;
            end
            q_objs = length(obj);
            score = zeros(q_objs,1);
            for j=1:q_objs
                q_repetitions = obj(j).repetitions;
    %             q_phot = cellfun(@length,[obj.T1_data.timetags(:)]);
                q_phot = zeros(q_repetitions,1);
                for i=1:q_repetitions
                    q_phot(i) = length(obj(j).T1_data.timetags{i,:});
                end

                % Get the mean value and standard deviation
                m_phot = mean(q_phot);
                s_phot = std(q_phot);

                % Set the thresholds
                a = 2; % Amount of standard deviations where to set the thresholds
                thr_low  = m_phot - a*s_phot;
                thr_high = m_phot + a*s_phot;

                % Find the runs bellow and beyond the thresholds
                run = 1:q_repetitions;
                filter_m = q_phot < thr_low | q_phot > thr_high;

                % Calculate the derivative of the run vector
                dq_phot = diff(q_phot);

                % Get statistics from the derivative
                mdq_phot = mean(dq_phot);
                sdq_phot = std(dq_phot);

                % Thresholds for the derivative
                b = a;
                dif_low  = mdq_phot - b*sdq_phot;
                dif_high = mdq_phot + b*sdq_phot;

                % Find the runs bellow and beyond the diff thresholds
                filter_dif = dq_phot < dif_low | dq_phot > dif_high;
                filter_dif = [filter_dif; 0] | [0; filter_dif];

                % force a range of data to remove
                range = [1:150];
                filter_f = zeros(q_repetitions,1);
                filter_f(range) = 1;

                % Integrate all the filtered data
                filter = filter_m | filter_dif | filter_f;
    %             filter = filter_dif;

                % Apply the filter. 1 -> remove;
                go_runs = run(~filter);
                nogo_runs = run(filter);

                % Save the filter in the object.
                obj(j).whitelist = go_runs;
                score(j) = length(go_runs) / q_repetitions;
                obj(j).QA_score = score(j);
                
                % Plot the data to filter
                if plt
                    figure
                    plot(go_runs, q_phot(go_runs),'.',nogo_runs, q_phot(nogo_runs),'.');
                    title(obj(j).name,'Interpreter','none');
                    xlabel('Run #')
                    ylabel('Quantity of Photons')
                end
            end
%             figure;
%             hist(q_phot,max(q_phot));
%             xlabel('Photons in a run')
%             ylabel('Quantity')

            %% Find edges in the data.
%             k = 10;
%             r = 40;
%             q_diff = length(dq_phot);
%             A = zeros(q_diff,1);
%             for i=1:(q_diff-k)
%                 A(i) = sum(abs(dq_phot(i:i+k)));
%             end
%             A = abs(A/k);
%             A = A - r;
%             jump = run(A>0);
%             hold on;
%             plot(run(jump),q_phot(jump),'+k');
%             figure;
%             plot(A,'.');
%             s = findchangepts(q_phot,'MaxNumChanges',34);
%             hold on;
%             plot(s,q_phot(s),'pk');
        end
      
      
        %% Reset the curve and the calculations.
        function reset_curve(obj)
            obj.T1_curve = []; % Integrated PL over the detection window of each reading pulse. Curve not normalised
            obj.T1_curve_n = []; % Normalised T1 curve.
            obj.T1_curve_f = []; % Fitted curve acording to the model used.
            obj.T1 = NaN; % Value of the spin - lattice relaxation constant.
            obj.PL_T1 = NaN; % Photoluminescence at T1.
            obj.fit_param = []; % Array with all the calculated parameters used to fit the model.
            obj.fit_error = NaN; % Norm of the residues of the fitting.
        end
        
        
        %% Plot the FFT of the integrated PL.
        function getFFT(obj,signal)
            % Check the optional parameter
            if exist('signal','var') == 0
                signal = obj.PL_int;
            end
            
            Ts = 10e-9; % Sampling period [s]
            Fs = 1/Ts; % Sampling frequency [Hz]
            L = length(signal);
            
            Y = fft(signal);
            P2 = abs(Y/L); %  two-sided spectrum
            P1 = P2(1:floor(L/2)+1); % single-sided spectrum
            P1(2:end-1) = 2*P1(2:end-1);
            
            f = Fs*(0:(L/2))/L; % Frequency axis [Hz]
            
            figure;
            plot(f,P1);
            title('Single-Sided Amplitude Spectrum of PL(t)');
            xlabel('f (Hz)');
            ylabel('|P1(f)|');
        end
        
        
        %% Print the name of all the files in the dataset
        function getname(obj)
           for i=1:length(obj)
            fprintf('%d\t%s\n',i,obj(i).name);
           end
        end
       
        
        %% Return the vector containing the darktimes used in the experiment.
        function dt = getdarktimes(obj)
            dt = obj.darktimes;
        end
        
        
        %% setWhitelist: sets the runs to be used to calculate the T1 curve.
        %                If the parameter 'list' is empty, then it resets the whitelist to 1:repetitions
        function setWhitelist(obj,list)
            if exist('list','var') == 0
                list = 1:obj.repetitions;
            end
            obj.whitelist = list;
        end
    end
end