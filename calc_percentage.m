clearvars -except T1_experiment;

% Remove repeated data
% I found in the water measurements that there is repeated files in the
% data.

nf = length(T1_experiment);
% order = [10 8 5 2 9 6 3 7 4 1]; % Water
% order = 1:10;
order = [10 8 5 2 9 6 3 7 4 1];

refe = 10;

%% ====== Mean T1 from all the T1s =======
% Find the existing groups
groups = {T1_experiment.group};
groups = intersect(groups, groups);
q_groups = length(groups);

mt1 = zeros(q_groups,1);
stdt1 = zeros(q_groups,1);
stdEt1 = zeros(q_groups,1);
rdst1 = zeros(q_groups,1);
q_sam = zeros(q_groups,1);
for gr=1:q_groups
    tmp = [];
    for i = 1:nf
        if ismember(T1_experiment(i).group,groups(gr))
           tmp = [tmp T1_experiment(i).T1];
        end
    end
    mt1(gr) = mean(tmp);
    stdt1(gr) = std(tmp); % Standard deviation
    q_sam(gr) = length(tmp);
    stdEt1(gr) = stdt1(gr)/sqrt(q_sam(gr)); % Standard error
    rdst1(gr) = stdt1(gr)/mt1(gr); % RSD    
end

Tb.Group = groups';
Tb.Tc2Exp = mt1; % us
Tb.RTc2Exp = 100 * mt1./ mt1(refe); % [%] water
% Tb.Tc2Exp = 100*mt1/mt1(10); % [%]
Tb.Std = stdt1; % us
Tb.StdErr = stdEt1; % us
Tb.RSD = 100 * rdst1; % [%]
Tb.RSE = 100 * stdEt1./mt1; % Relative standard error
Tb.Samples = q_sam;
Tb = struct2table(Tb);
disp(Tb);
err_neg = Tb.RSE;
err_pos = err_neg;



% Reference all the T1s to the mean value of the water samples
for gr=1:q_groups
    tmp = [];
    for i = 1:nf
        if ismember(T1_experiment(i).group,groups(gr))
           tmp = [tmp 100*T1_experiment(i).T1/mt1(refe)]; % Referencing to water
        end
    end
    Rmt1(gr) = mean(tmp);
    Rstdt1(gr) = std(tmp); % Standard deviation
    q_sam(gr) = length(tmp);
    RstdEt1(gr) = Rstdt1(gr)/sqrt(q_sam(gr)); % Standard error
    Rrdst1(gr) = Rstdt1(gr)/Rmt1(gr); % RSD    
end

RTb.Group = groups';
RTb.RTc2Exp = Rmt1'; % [%]
RTb.Std = Rstdt1'; % [%]
RTb.StdErr = RstdEt1'; % [%]
RTb.Samples = q_sam;
RTb = struct2table(RTb);
disp(RTb);
err_neg = RTb.StdErr;
err_pos = err_neg;


% return

% Create the figure
figA = figure;
setFigureColors(figA);

conc = 1:10;% nM

errorbar(conc, RTb.RTc2Exp(order), err_neg(order), err_pos(order),'.');

title('Water - % water, SE');
xlabel('Concentration of Gd^{3+} [nM]');
ylabel('T_1 %');
xticks(conc);
xticklabels({'0' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'});
grid on;
formatFigure(figA);
% 
% saveFigure(figA,'fig4_Water');


% plot(Tb.Tc2Exp(order),'.');

return
%% Fill the figure
errorbar(conc, corona_groups.T1s(concCO), corona_groups.err_neg(concCO), corona_groups.err_pos(concCO),'.');
hold on;
errorbar(conc, water_groups.T1s(concWG), water_groups.err_neg(concWG), water_groups.err_pos(concWG),'.');
% errorbar(conc, PBSGroups.T1s(concPG), PBSGroups.err_neg(concWG), PBSGroups.err_pos(concWG),'.');
errorbar(conc(concPG), PBS_groups.T1s(concPG), PBS_groups.err_neg(concWG(concPG)), PBS_groups.err_pos(concWG(concPG)),'.');

xlabel('Concentration of Gd^{3+} [nM]');
ylabel('T_1 %');
ylim([0 150]);
xlim([0.5 10.5]);
xticks(conc);
xticklabels({'0' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'});
legend({'Corona','Water','PBS'},'Location','best','Orientation','horizontal');
% figA.Children(2).XScale = 'log';
grid on;
formatFigure(figA);
% 
saveFigure(figA,'fig4_Water-PBS-FBS_RSD');
