function [pics_t, pics_f] = pics_spectraux(S, eta_t, eta_f, epsilon)
    se = ones(eta_t, eta_f);
    S_max = imdilate(S, se);
    pos_eg = S_max == S;
    pos_eps = S_max >= epsilon;
    pos = pos_eg .* pos_eps;
    [pics_f, pics_t, ~] = find(pos);
end