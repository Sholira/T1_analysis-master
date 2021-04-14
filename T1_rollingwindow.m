function [segment] = T1_rollingwindow(T1_experiment, filepath, model, window, shift, time, NORM, combi, bootsam)
nf = length(T1_experiment);
t = linspace(0,time,nf);
bootcom = round(bootsam/4);


for j = 1:nf
    
    segment = floor(1 + (T1_experiment(j).repetitions-window)/shift);
    fromto(:,1) = 1:shift:(T1_experiment(j).T1_data.repetitions - window + 1);
    fromto(:,2) = window:shift:T1_experiment(j).T1_data.repetitions;
    t = linspace(0,time,segment);

    timetags = T1_experiment(j).T1_data.timetags;
    triggers = T1_experiment(j).T1_data.triggers;
    corr = T1_experiment(j).correction;
    r_width = T1_experiment(j).read_width;
    darktimes = T1_experiment(j).darktimes;
    curve = zeros(segment,length(T1_experiment(j).darktimes));
    parfor x = 1:segment
        seg_timetags = timetags(fromto(x,1):fromto(x,2));
        seg_triggers = triggers(fromto(x,1):fromto(x,2));
        tmp = getT1curve(seg_triggers, seg_timetags, corr, r_width);
        if NORM == true
            tmp = doNORM(tmp,'tail');
        end
        curve(x,:) = tmp;
        tmp_sr = doNORM(tmp,'stretch');
        curve_sr(x,:) = tmp_sr;
    end
    
    if combi == true
        bootcurvewin = strcat(filepath, T1_experiment(j).name,'_RW_fol\');
        for x = 1:segment
            seg_timetags = timetags(fromto(x,1):fromto(x,2));
            seg_triggers = triggers(fromto(x,1):fromto(x,2));
            [mod_biL(x,:), mod_biW(x,:), mod_si(x,:), mod_sr(x,:)] = T1_bootstrap(strcat(T1_experiment(j).name,'_window_',num2str(x)),seg_timetags, seg_triggers, darktimes,  T1_experiment(j).read_width, T1_experiment(j).correction, model, NORM, bootcom, bootcurvewin);
        end
    end
    
    parnam = T1_experiment(j).name;
    
    if strcmp(model, 'biexp') || strcmp(model, 'all')
        T1curve_bi = zeros(size(curve));
        fitparam_bi = zeros(segment, 5);
        T1L = zeros(segment,1);
        T1W = zeros(segment,1);
        PLT1_bi = zeros(segment,1);
        fiterror_bi = zeros(segment,1);
        for x = 1:segment
            if x == 1
                [T1curve_bi(x,:), fitparam_bi(x,:), PLT1_bi(x), fiterror_bi(x)] = biexp_model(darktimes, curve(x,:));
                [T1L(x,1), T1W(x,1)] = select_param(fitparam_bi(x,:));
            else
                [T1curve_bi(x,:), fitparam_bi(x,:), PLT1_bi(x), fiterror_bi(x)] = biexp_model(darktimes, curve(x,:), fitparam_bi((x-1),:));
                [T1L(x,1), T1W(x,1)] = select_param(fitparam_bi(x,:));
            end
        end
        docnam = strcat(filepath,parnam,'_RW_biexp.xlsx');
        T = table(t', T1L, T1W, fitparam_bi, 'VariableNames', {'time', 'T1_long', 'T1_weight', 'parameters'});
        writetable(T,docnam)
        clear T docnam
        if combi == true
            docnam2 = strcat(filepath,parnam,'_RW_biexp_boot.xlsx');
            T = table(t', mod_biL, mod_biW, 'VariableNames', {'time', 'T1_mode_long', 'T1_mode_weight'});
            writetable(T,docnam2);
            clear T docnam2
        end
    end
    
    if strcmp(model,'siexp') || strcmp(model,'all')
        T1curve_si = zeros(size(curve));
        fitparam_si = zeros(segment, 3);
        PLT1_si = zeros(segment, 1);
        fiterror_si = zeros(segment, 1);
        for x = 1:segment
            if x ==1
                [T1curve_si(x,:), fitparam_si(x,:), PLT1_si(x), fiterror_si(x)] = singexp_model(darktimes, curve(x,:));
            else
                [T1curve_si(x,:), fitparam_si(x,:), PLT1_si(x), fiterror_si(x)] = singexp_model(darktimes, curve(x,:), fitparam_si((x-1),:));
            end
        end
        T1S = fitparam_si(:,3);
        docnam = strcat(filepath,parnam,'_RW_siexp.xlsx');
        T = table(t', T1S, fitparam_si, 'VariableNames', {'time', 'T1_value', 'parameters'});
        writetable(T,docnam)
        clear T docnam
        if combi == true
            docnam2 = strcat(filepath, parnam, '_RW_siexp_boot.xlsx');
            T = table(t', mod_si, 'VariableNames', {'time', 'T1_mode'});
            writetable (T, docnam2)
            clear T docnam2
        end
    end
    
    if strcmp(model,'srexp') || strcmp(model,'all')
        T1curve_st = zeros(size(curve));
        T1R = zeros(segment, 1);
        PLT1_st = zeros(segment, 1);
        fiterror_st = zeros(segment, 1);
        for x = 1:segment
            if x == 1
                [T1curve_st(x,:), T1R(x,:), PLT1_st(x), fiterror_st(x)] = stretch_model(darktimes, curve_sr(x,:));
            else
                [T1curve_st(x,:), T1R(x,:), PLT1_st(x), fiterror_st(x)] = stretch_model(darktimes, curve_sr(x,:), T1R((x-1),:));
            end
        end
        docnam = strcat(filepath, parnam, '_RW_stretch.xlsx');
        T = table (t', T1R, 'VariableNames', {'time', 'T1_value'});
        writetable(T,docnam)
        clear T
        if combi == true
            docnam2 = strcat(filepath, parnam, '_RW_stretchexp_boot.xlsx');
            T = table(t', mod_sr, 'VariableNames', {'time', 'T1_mode'});
            writetable(T, docnam2)
            clear T docnam2
        end
    end
                
end

end