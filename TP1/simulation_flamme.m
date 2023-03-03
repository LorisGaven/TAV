clear;
close all;

taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

load donnees;
beta_0 = bord_inf(1,1);
gamma_0 = bord_sup(1,1);

load parametres_lois;

load texture;
[nb_lignes,nb_colonnes] = size(texture);

% Simulation d'une flamme de bougie :
figure('Name','Simulation d''une flamme de bougie','Position',[0.4*L,0.05*H,0.6*L,0.7*H]);
l_image = 1000;					% Largeur de l'image
h_image = 1000;					% Hauteur de l'image
h_flamme = round(0.85*h_image);			% Hauteur de la flamme
x = transpose(0:1/(h_flamme-1):1);		% Abscisses normalisees entre 0 et 1
nb_images = 100;				% Longueur de la sequence
for j = 1:nb_images
	I = zeros(h_image,l_image);
	[y_inf,y_sup] = tirage_aleatoire(x,moyennes,ecarts_types,beta_0,gamma_0);

	if sum(y_inf<y_sup)==0
		for ligne = 1:h_flamme
			ligne_tex = round((nb_lignes*(h_flamme-ligne)+ligne-1)/(h_flamme-1));
			colonne_min = floor(l_image/2+(y_sup(ligne)-(beta_0+gamma_0)/2)*l_image/150);
			colonne_max = ceil(l_image/2+(y_inf(ligne)-(beta_0+gamma_0)/2)*l_image/150);
			largeur_flamme = colonne_max-colonne_min;
			for colonne = max(colonne_min,1):min(colonne_max,l_image)
				colonne_tex = round((colonne-colonne_min)*(nb_colonnes-1)/largeur_flamme+1);
				I(ligne,colonne) = round(255*texture(ligne_tex,colonne_tex));
			end
		end
		imagesc(I);
		axis xy;
		axis off;
		colormap(hot);

		pause(0.1);
	end
end

function [moyennes,ecarts_types] = estimation_lois(liste_parametres)
    nb_param = size(liste_parametres,2);
    moyennes = sum(liste_parametres,2) ./ nb_param;
    ecarts_types = sqrt(sum((liste_parametres - moyennes).^2,2) ./ nb_param);
end

function b = bernstein(d,k,x)
    b = nchoosek(d,k) * x^k * (1-x)^(d-k);
end

function y = bezier(x, param_0, param)
    p = length(x);
    d = length(param);
    y = zeros(p,1);
    for i = 1:p
        y(i) = param_0 * bernstein(d,0,x(i));
        for j = 1:d
            y(i) = y(i) +  param(j) * bernstein(d, j, x(i));
        end
    end
end

function [y_inf,y_sup] = tirage_aleatoire(x,moyennes,ecarts_types,beta_0,gamma_0)
   echant = moyennes + ecarts_types .* randn(length(moyennes),1);
   d = 5;
   betas = [echant(1:(d-1)); echant(2*d-1)];
   gammas = echant(d:(2*d-1));
   y_sup = bezier(x, gamma_0, gammas);
   y_inf = bezier(x, beta_0, betas);
end
