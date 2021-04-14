function curve = getT1curve(triggers,timetags,correction,read_width)
    % Calculates the start and stop flags of the detection window.
    [starts, stops] = find_start_stop(triggers,correction,read_width);

    % Count the photons inside the detection window.
    detections = cellfun(@count_photons,timetags,starts,stops,'UniformOutput', false);
    detections = cell2mat(detections);

    curve = sum(detections);

    % Remove the detection of the first polarizing pulse.
    curve = curve(2:end);
end