function y_modifie = etirement_temporel(y, f_ech, accel)
    Y = TFCT(y, f_ech, 2048, 512, "hann");
    C = 1:accel:size(Y,2);
    phi = angle(Y(:,1));
    y_modifie = zeros(size(Y,1), size(C,2));
    Y = padarray(Y, [0,1], 0, 'post');
    for i = 1:length(C)
        c = floor(C(i));
        alpha = C(i) - c;
        rho = (1 - alpha) .* abs(Y(:,c)) + alpha .* abs(Y(:,c+1));
        y_modifie(:,i) = rho .* exp(j*phi);
        d_phi = angle(Y(:,c+1)) - angle(Y(:,c));
        phi = phi + d_phi;
    end
    y_modifie = ITFCT(y_modifie, f_ech, 512, "hann");
end