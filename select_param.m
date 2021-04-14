function [T1L, T1W] = select_param(fitparam)

T1L = max([fitparam(3) fitparam(5)]);
[~, ind] = max([fitparam(2) fitparam(4)]);
if ind == 1
    T1W = fitparam(3);
else
    T1W = fitparam(5);
end
end
