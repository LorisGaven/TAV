function paires = appariement(pics_t, pics_f, n_v, delta_t, delta_f)
    n = size(pics_t, 1);
    paires = [];
    for i = 1:n
        nb_paires = 0;
        for j = 1:n
            m_diff = pics_t(j) - pics_t(i);
            k_diff = abs(pics_f(i) - pics_f(j));
            if nb_paires >= n_v
                break;
            end
            if m_diff > 0 && m_diff <= delta_t && k_diff <= delta_f && i ~= j
                paires = [paires; pics_f(i) pics_f(j) pics_t(i) pics_t(j)];
                nb_paires = nb_paires + 1;
            end
        end
    end
end