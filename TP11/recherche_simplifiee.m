function resultats = recherche_simplifiee(identifiants, bdd)
    resultats = [];
    for i = 1:length(identifiants)
        if bdd.isKey(identifiants(i))
            entrees = bdd(identifiants(i));
            resultats = [resultats; entrees(:,2)];
        end
    end
end