function [starts, stops] = find_start_stop(triggers, correction, read_width)
    % Convert correction and read_with [us] to ten of nanoseconds (one conunt of the DAQ
    % counter takes 10ns).
    correction = correction * 100;
    read_width = read_width * 100;
    
    % Add a first trigger at time = 0. This corresponds to the first
    % polarizating pulse.
    triggers = cellfun(@(x)([0 x]),triggers,'UniformOutput', false);

    % Shift every trigger by the correction.
    triggers = cellfun(@(x)(x + correction),triggers,'UniformOutput', false);
    
    starts = triggers;
    
    % Calculates the sop flags. 
    stops = cellfun(@(x)(x + read_width),triggers,'UniformOutput', false);
end