function detections = count_photons(timetags,starts,stops)
    q_pulses = length(starts);
    detections = zeros(1, q_pulses);
    for i=1:q_pulses
        wind = timetags(timetags>=starts(i) & timetags<=stops(i));
        detections(i) = length(wind);
    end
end