function [T1s] = T1_segmentation(T1_experiment, segment, model, filepath, NORM)
nf = length(T1_experiment);

for j = 1:nf
    q_run_seg = floor(T1_experiment(j).T1_data.repetitions/segment);
    remanent = T1_experiment(j).T1_data.repetitions - q_run_seg*segment;
    if remanent == 0
        q_run_seg = q_run_seg * ones(1,segment);
    else
        q_run_seg = [q_run_seg * ones(1,segment-1), q_run_seg + remanent];
    end
    
    seg_start = zeros(1,segment);
    seg_stop = zeros(1,segment);
    seg_stop = 0;
    for x = 1:segment
        seg_start(x) = seg_stop + 1;
        seg_stop_tmp = sum(q_run_seg(1:x));
        seg_stop(1,x) = seg_stop_tmp;
    end
    fromto(:,1) = seg_start;
    fromto(:,2) = seg_stop;
    
    curve = zeros(segment, length(T1_experiment(j).darktimes));
    for x = 1:segment
        seg_timetags = T1_experiment(j).T1_data.timetags(seg_start(x):seg_stop(x));
        seg_triggers = T1_experiment(j).T1_data.triggers(seg_start(x):seg_stop(x));
        tmp = getT1curve(seg_triggers, seg_timetags, T1_experiment(j).correction, T1_experiment(j).read_width);
        if NORM == true
            tmp = doNORM(tmp,'tail');
        end
        curve(x,:) = tmp;
        tmp_sr = doNORM(tmp,'stretch');
        curve_sr = tmp_sr;
    end
    
    parnam = T1_experiment(j).name;
    if strcmp(model,'biexp') || strcmp(model,'all')
        T1curve_bi = zeros(size(curve));
        fitparam_bi = zeros(segment, 5);
        T1L = zeros(segment,1);
        T1W = zeros(segment,1);
        PLT1_bi = zeros(segment, 1);
        fiterror_bi = zeros(segment,1);
        for x = 1:segment
            [T1curve_bi(x,:), fitparam_bi(x,:), PLT1_bi(x), fiterror_bi(x)] = biexp_model(T1_experiment(j).darktimes,curve(x,:));
            T1L(x,1) = max([fitparam_bi(x,3) fitparam_bi(x,5)]);
            [~, ind] = max([fitparam_bi(x,2) fitparam_bi(x,4)]);
            if ind == 1
                T1W(x,1) = fitparam_bi(x,3);
            else
                T1W(x,1) = fitparam_bi(x,5);
            end
        end
        docnam = strcat(filepath,parnam,'_segment_biexp.xlsx');
        T = table (T1L, T1W, fitparam_bi, 'VariableNames', {'T1_long', 'T1_weight', 'parameters'});
        writetable(T,docnam)
        clear T
    end
    
    if strcmp(model,'siexp') || strcmp(model,'all')
        T1curve_si = zeros(size(curve));
        fitparam_si = zeros(segment, 3);
        T1S = zeros(segment, 1);
        PLT1_si = zeros(segment, 1);
        fiterror_si = zeros(segment, 1);
        for x = 1:segment
            [T1curve_si(x,:), fitparam_si(x,:), PLT1_si(x), fiterror_si(x)] = singexp_model(T1_experiment(j).darktimes,curve(x,:));
            T1S(x,1) = fitparam_si(x,3);
        end
        docnam = strcat(filepath,parnam,'_segment_siexp.xlsx');
        T = table (T1S, fitparam_si, 'VariableNames', {'T1_value', 'parameters'});
        writetable(T,docnam)
        clear T
    end
    
    if strcmp(model,'srexp') || strcmp(model,'all')
        T1curve_st = zeros(size(curve));
        fitparam_st = zeros(segment, 1);
        T1R = zeros(segment, 1);
        PLT1_st = zeros(segment, 1);
        fiterror_st = zeros(segment, 1);
        for x = 1:segment
            [T1curve_st(x,:), fitparam_st(x,:), PLT1_st(x), fiterror_st(x)] = stretch_model(T1_experiment(j).darktimes,curve_sr(x,:));
            T1R = fitparam_st(x,1);
        end
        docnam = strcat(filepath,parnam,'_segment_stretch.xlsx');
        T = table (T1R, 'VariableNames', {'T1_value'});
        writetable(T,docnam)
        clear T
    end
            
        
    

    clear T1curve_f PLT1 fiterror curve
end
end