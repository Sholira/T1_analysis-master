function value = get_value(text_line)
% Extract one numerical value from a string. It is assumed a
% name=value format. So the number is read from the string after the
% '=' sign.
%
    pos = strfind(text_line,'=') + 1; % find the position after the '=' sign.
    tmp = text_line(pos:end); % get the number as a string.
    value = sscanf(tmp,'%f',1); % convert the string to a Real number.
end