function [F_h, F_p] = HPSS(Y)
    F_h = medfilt2(Y, [1, 17]);
    F_p = medfilt2(Y, [17, 1]);
end