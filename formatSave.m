function formatSave(fig)
    % Formating the figure for saving
    % Parameters
    file_width = 14;
    file_height = 8.6;
    alw = 2; % AxesLineWidth
    fsz = 12; % Fontsize
    plw = 1.2; % Plot line width
    pmz = 10; % Plot Marker size

    % Configure
    ax = gca;
    ax.FontSize = fsz;
    % ax.LineWidth = alw;
    ax.Box = 'off';
    ax.FontWeight = 'bold';
    set(findall(gca, 'Type', 'Line'),'LineWidth',plw);
    set(findall(gca, 'Type', 'Line'),'MarkerSize',pmz);

    fig.Units = 'centimeters';
    fig.Position = [0 0 file_width file_height];
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 file_width file_height];

%     legend('off');
%     legend({'Gadolinium','Water'});

    % === Save figure ===
    % Save the image in the same folder where the data files are
    f_name = [file_path f_name '.png'];
    
    % Save relaxation curves
%     print(figA,'-dpng','-r300',f_name);

    % Save time constants plot
    print(figB,'-dpng','-r300',f_name);
    
%     f_name = sprintf('rational_function');
%     f_name = sprintf('comparation');
%     print(figC,'-dpng','-r300',f_name);   
end