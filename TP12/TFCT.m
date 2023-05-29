function [Y, valeurs_t, valeurs_f] = TFCT(y, f_ech, n_fenetre, n_decalage, fenetre)

    L = length(y);
    
    if strcmp(fenetre, 'hann')
        w = hann(n_fenetre);
    elseif strcmp(fenetre, 'rect')
	    w = ones(n_fenetre, 1);
    end
    
    y = buffer(y, n_fenetre,n_fenetre - n_decalage,"nodelay");
    
    tfct = fft(y.*w);
    
    Y = tfct(1:(n_fenetre / 2) + 1,:);
    
    k = 0:(n_fenetre/2);
    m = 0:(floor((L - n_fenetre) / n_decalage) + 1);
    valeurs_t = m*n_decalage/ f_ech;
    valeurs_f = k*f_ech/n_fenetre;
end