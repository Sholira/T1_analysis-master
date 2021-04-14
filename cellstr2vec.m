function arr_num = cellstr2vec(cellsrt)
    arr_num = str2double(strsplit(cellsrt));
    arr_num = arr_num(arr_num~=0);
end