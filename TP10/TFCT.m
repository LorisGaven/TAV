function [Y, valeurs_t, valeurs_f] = TFCT(y, f_ech, n_fenetre, n_decalage, fenetre)
    N = n_fenetre;
    L = length(y);
    w = hann(n_fenetre);
    buf = buffer(y, n_fenetre, n_fenetre-n_decalage, "nodelay");
    if fenetre == "hann"
        yw = buf .* w;
    else
        yw = buf;
    end
    Y = fft(yw);
    Y = Y(1:(floor(n_fenetre/2)+1),:);
    k = 0:(floor(N/2));
    m = 0:(floor((L-N)/n_decalage)+1);
    valeurs_t = (m * n_decalage) / f_ech;
    valeurs_f = (k * f_ech) / N;
end