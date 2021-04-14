function T1_data_norm = doNORM(T1_curve,type)
    % There are three types of normalisation:
    % 'tail' : The last values of the data are normalised to 1.
    % 'head' : The first values of the data are normalised to 1.
    % 'stretch' : The last values are normalised to 0 and the first values
    %             to 1.

    % Parameters
    tail_size = 10;
    head_size = 3;

    % Normalisation
    if strcmp(type,'stretch')
        % Find the mean at the ends of the curve.
        mean_tail = mean(T1_curve(end-tail_size:end));            

        % Move the curve down, setting the tail at zero
        T1_curve = (T1_curve - mean_tail);

        % Find the mean at the begining of the curve.
        mean_head = mean(T1_curve(1:head_size));

        % Scale the curve to the reange [0, 1]
        T1_curve = T1_curve / mean_head;

    elseif strcmp(type,'head')
        % Find the mean at the begining of the curve.
        mean_head = mean(T1_curve(1:head_size));

        % Move the curve, setting the head to one
        T1_curve = T1_curve / mean_head;

    % By default, normalise the tail of the curve at one.
    else
        % Find the mean at the ends of the curve.
        mean_tail = mean(T1_curve(end-tail_size:end)); 

        % Move the curve down, setting the tail to one
        T1_curve = T1_curve / mean_tail;
    end

%             % Normalyse the standard deviation
%             if isfield(T1_data(i),'PL_std')
%                 T1_data(i).PL_std = T1_data(i).PL_std / mean_tail;
%             end
%         end

% Set the return value
T1_data_norm = T1_curve;
end
