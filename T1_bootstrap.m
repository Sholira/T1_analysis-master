function [mod_biL, mod_biW, mod_si, mod_sr] = T1_bootstrap(name, timetags, triggers, darktimes, rwid, corr, model, NORM, bootsam, filepath)
%% this function creates a bootstrap for the timetags and triggers that are applied to the system.
% To speed up the process the program writes away a lot of .mat files, so
% when running this script enough memory has to be free.
% The high number of inputs allows to script to be used more freely and
% omitting data that is not needed for the bootstrap itself.
    
%% make the folder for the bootstrap samples

bootnam = strcat(filepath, name,'_bootsam\');
if isfolder(bootnam) ==0
    mkdir (bootnam)
end

%%  Create the bootsrap samples and calculate the T1

rep = length(timetags);

bootcurvfol = strcat(bootnam,name, '_bootcurve\');
if isfolder(bootcurvfol) ==0
    mkdir(bootcurvfol);
end


parfor j = 1:bootsam
    rng('shuffle');
    bootsub = randi(rep,[1 rep]);
    timet = timetags(bootsub);
    trigg = triggers(bootsub);
    bootsub = bootsub';
    tmp = getT1curve(trigg,timet,corr,rwid);
    boot = strcat(bootnam,num2str(j),'.mat');
    parsave(boot,bootsub);
    
    if strcmp(model,'biexp') || strcmp(model,'all')
        if NORM == true
            tmp = doNORM(tmp,'tail');
        end
        [T1curve_bi, fitparam_bi(j,:), ~, ~] = biexp_model(darktimes,tmp);
        bootcurve = strcat(bootcurvfol,name,'_biexp_curve_',num2str(j),'.mat');
        parsave(bootcurve,T1curve_bi);
    end
    
    if strcmp(model,'siexp') || strcmp(model,'all')
        if NORM == true
            tmp = doNORM(tmp,'tail')
        end
        [T1curve_si, fitparam_si(j,:), ~, ~] = singexp_model(darktimes,tmp);
        bootcurve = strcat(bootcurvfol,name,'_siexp_curve_',num2str(j),'.mat');
        parsave(bootcurve, T1curve_si);
    end
    
    if strcmp(model,'srexp') || strcmp(model,'all')
        tmp = doNORM(tmp,'stretch');
        [T1curve_sr, fitparam_sr(j,:), ~, ~] = stretch_model(darktimes,tmp);
        bootcurve = strcat(bootcurvfol,name,'_srexp_curve_',num2str(j),'.mat');
        parsave(bootcurve, T1curve_sr);
    end
end

if strcmp(model,'biexp') || strcmp(model,'all')
    T1L = fitparam_bi(:,5);
    for k=1:length(T1L)
        if fitparam_bi(k,2) > fitparam_bi(k,4)
            T1W(k,1) = fitparam_bi(k,3);
        else
            T1W(k,1) = fitparam_bi(k,5);
        end
    end
    T1_biL = strcat(filepath,name,'_boot_biexp-long_T1-list.xlsx');
    T = table(T1L, 'VariableNames',{'T1_long'});
    writetable(T, T1_biL)
    T1_biW = strcat(filepath,name,'_boot_biexp-weight_T1-list.xlsx');
    Y = table(T1W, 'VariableNames',{'T1_weight'});
    writetable(Y,T1_biW)
    clear T1_biL T T1_biW Y
end

if strcmp(model,'siexp') || strcmp(model,'all')
    T1S = fitparam_si(:,3);
    T1_si = strcat(filepath,name, '_boot_siexp_T1-list.xlsx');
    S = table(T1S, 'VariableNames', {'T1_siexp'});
    writetable(S, T1_si)
    clear S T1_si
end

if strcmp(model,'srexp') || strcmp(model,'all')
    T1R = fitparam_sr;
    T1_sr = strcat(filepath, name, '_boot_srexp_T1-list.xlsx');
    R = table(T1R, 'VariableNames', {'T1_srexp'});
    writetable(R, T1_sr)
    clear R T1_sr
end

clear bootsub tmp bootcurve T1curve_bi T1curve_si T1curve_sr trigg timet

%% create the T1 distribution based on KDE

min95 = 0.025;
max95 = 0.975;

if strcmp(model,'biexp') || strcmp(model,'all')
    maxhis_biL = max(T1L);
    meanhis_biL = mean(T1L);
    if maxhis_biL > (meanhis_biL*5)
        maxhis_biL = meanhis_biL*5;
    end
    [~, density, xmesh, cdf] = kde (T1L, 2^12, 0, maxhis_biL);
    x = find(density == max(density), 1, 'first');
    maxx_bi = xmesh(x);
    Y = [xmesh' density];
    writetable(table(Y), char(strcat(filepath,name,'_biexp-long_kde.xlsx')));
    min_biL = find(cdf >= min95, 1, 'first');
    min95_biL = xmesh(min_biL);
    max_biL = find(cdf >= max95, 1, 'first');
    max95_biL = xmesh(max_biL);
    mod_biL = [maxx_bi min95_biL max95_biL];
    clear cdf
    
    maxhis_biW = max(T1W);
    meanhis_biW = mean(T1W);
    if maxhis_biW > (meanhis_biW*5)
        maxhis_biW = meanhis_biW*5;
    end
    [~, density, xmesh, cdf] = kde (T1W, 2^12 ,0 , maxhis_biW);
    y = find(density == max(density), 1, 'first');
    maxx_biw = xmesh(y);
    Z = [xmesh' density];
    writetable(table(Z), char(strcat(filepath, name, '_biexp-weight_kde.xlsx')));
    min_biW = find(cdf >= min95, 1, 'first');
    min95_biW = xmesh(min_biW);
    max_biW = find(cdf >= max95, 1, 'first');
    max95_biW = xmesh(max_biW);
    mod_biW = [maxx_biw min95_biW max95_biW];
    clear maxhis_biL maxhis_biW meanhis_biL meanhis_biW x maxx_bi U Z y maxx_biw Z Q density xmesh cdf
end

if strcmp(model,'siexp') || strcmp(model,'all')
    maxhis_si = max(T1S);
    meanhis_si = mean(T1S);
    if maxhis_si > (meanhis_si*5)
        maxhis_si = meanhis_si*5;
    end
    [~, density, xmesh, cdf] = kde(T1S, 2^12, 0, maxhis_si);
    x = find(density == max(density), 1, 'first');
    maxx_si = xmesh(x);
    U = [xmesh' density];
    writetable(table(U), char(strcat(filepath, name, '_siexp_kde.xlsx')))
    min_si = find(cdf >= min95, 1, 'first');
    min95_si = xmesh(min_si);
    max_si = find(cdf >= max95, 1, 'first');
    max95_si = xmesh(max_si);
    mod_si = [maxx_si min95_si max95_si];
    clear maxhis_si meanhis_si x density xmesh  maxx_si U P cdf
end

if strcmp(model, 'srexp') || strcmp(model, 'all')
    maxhis_sr = max(T1R);
    meanhis_sr = mean(T1R);
    if maxhis_sr > (meanhis_sr*5)
        maxhis_sr = meanhis_sr*5;
    end
    [~, density, xmesh, cdf] = kde(T1R, 2^12, 0, maxhis_sr);
    x = find(density == max(density), 1, 'first');
    maxx_sr = xmesh(x);
    A = [xmesh', density];
    writetable(table(A), char(strcat(filepath, name, '_srexp_kde.xlsx')))
    min_sr = find(cdf >= min95, 1, 'first');
    min95_sr = xmesh(min_sr);
    max_sr = find(cdf >= max95, 1, 'first');
    max95_sr = xmesh(max_sr);
    mod_sr = [maxx_sr min95_sr max95_sr];
    clear maxhis_sr meanhis_sr density xmesh x maxx_sr A S cdf
end

    

%clear T bootfile bootsub timetags triggers bootsub tmp bootfile T boot

end