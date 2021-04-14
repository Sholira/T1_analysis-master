function [ ] = model_fit(T1_experiment, filepath, model, NORM)
% this script takes the data and plots the data based on the input
clear y_fitted params y_TC resnorm
nf = length(T1_experiment);

if strcmp(model, 'biexp') || strcmp(model,'all')
    for j = 1:nf
        dt = T1_experiment(j).darktimes;
        if NORM == true
            y1 = T1_experiment(j).T1_curve_n;
            y = doNORM(y1,'tail');
        else
            y = T1_experiment(j).T1_curve_n;
        end
        [y_fitted(j,:), params(j,:), y_TC(j), resnorm(j)] = biexp_model(dt, y);
        T1L(j,1) = max([params(j,3) params(j,5)]);
        [~, ind] = max([params(j,2) params(j,4)]);
        if ind == 1
            T1W(j,1) = params(j,3);
        else
            T1W(j,1) = params(j,5);
        end
        parnam{j,1} = T1_experiment(j).name;
    end
    docnam = strcat(filepath, 'T1_list_biexponential.xlsx');
    T = table(parnam, T1L, T1W, params, y_TC', resnorm', 'VariableNames', {'filenam', 'T1_long', 'T1_weight', 'parameters', 'contrast', 'resnorm'});
    writetable(T,docnam)
    clear T y_fitted params y_TC resnorm j
end
    
    
if strcmp(model, 'siexp') || strcmp(model,'all')
    for j = 1:nf
        dt = T1_experiment(j).darktimes;
        if NORM == true
            y1 = T1_experiment(j).T1_curve_n;
            y = doNORM(y1,'head');
        else
            y = T1_experiment(j).T1_curve_n;
        end
        [y_fitted(j,:), params(j,:), y_TC(j), resnorm(j)] = singexp_model(dt,y);
        parnam{j,1} = T1_experiment(j).name;
    end
    docnam = strcat(filepath, 'T1_list_singleexponential.xlsx');
    T = table(parnam, params, y_TC', resnorm', 'VariableNames', {'filenam', 'parameters', 'contrast', 'resnorm'});
    writetable(T,docnam)
    clear T y_fitted params y_TC resnorm j 
end
    
    
if strcmp(model, 'srexp') || strcmp(model, 'all')
    for j = 1:nf
        dt = T1_experiment(j).darktimes;
        y1 = T1_experiment(j).T1_curve_n;
        y = doNORM(y1, 'stretch');
        [y_fitted(j,:), params(j,:), y_TC(j), resnorm(j)] = stretch_model(dt,y);
        parnam{j,1} = T1_experiment(j).name;
    end
    docnam = strcat(filepath, 'T1_list_stretchexponential.xlsx');
    T = table(parnam, params, y_TC', resnorm', 'VariableNames', {'filenam', 'parameters', 'contrast', 'resnorm'});
    writetable(T,docnam)
    clear T y_fitted params y_TC resnorm j
end


end