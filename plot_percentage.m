% load('X:\Data\water.mat')
% load('X:\Data\PBS.mat')
% load('X:\Data\FBS.mat')

% Create the figure
figA = figure;
setFigureColors(figA);

% Concetrations
% conc = [0.1 1 10 100 1e3 10e3 100e3 1e6 10e6 100e6];% nM
conc = 1:10;% nM
concWG = [10 8 5 2 9 6 3 7 4 1];
concPG = 1:5;
concCO = [10 8 5 2 9 6 3 7 4 1];

%% Fill the figure
errorbar(conc, RTb_water.RTc2Exp(concWG), RTb_water.StdErr(concWG), RTb_water.StdErr(concWG),'.');
hold on;
errorbar(conc, RTb_FBS.RTc2Exp(concCO), RTb_FBS.StdErr(concCO), RTb_FBS.StdErr(concCO),'.');
errorbar(conc(concPG), RTb_PBS.RTc2Exp(concPG), RTb_PBS.StdErr(concPG), RTb_PBS.StdErr(concPG),'.');


xlabel('Concentration of Gd^{3+} [nM]');
ylabel('T_1 %');
ylim([0 150]);
xlim([0.5 10.5]);
xticks(conc);
xticklabels({'0' '10^0' '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'});
legend({'Water','Corona','PBS'},'Location','best','Orientation','horizontal');
grid on;
formatFigure(figA);

saveFigure(figA,'fig4_Water-PBS-FBS_perc');