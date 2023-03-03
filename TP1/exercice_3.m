clear;
close all;

taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

load donnees;
beta_0 = bord_inf(1,1);
gamma_0 = bord_sup(1,1);

load points_de_controle;

% Estimation des lois normales :
[moyennes,ecarts_types] = estimation_lois(liste_parametres);

% Simulation de silhouettes par tirages aleatoires :
figure('Name','Tirage aleatoire de silhouettes','Position',[0.4*L,0.05*H,0.6*L,0.7*H]);
p = size(bord_inf,1);
x = transpose(0:1/(p-1):1);
nb_images = 50;				% Longueur de la sequence simulee
for j = 1:nb_images
	[y_inf,y_sup] = tirage_aleatoire(x,moyennes,ecarts_types,beta_0,gamma_0);

	if sum(y_inf<y_sup)==0
		plot(x,y_inf,'Color','r','LineWidth',2);
		hold on;
		plot(x,y_sup,'Color','r','LineWidth',2);
		axis([0,1.01,60,150]);
		axis ij;
		set(gca,'FontSize',20);
		xlabel('$x$','FontSize',30,'Interpreter','Latex');
		ylabel('$y$','FontSize',30,'Interpreter','Latex','Rotation',0);

		pause(0.2);
		hold off;
	end
end

save parametres_lois moyennes ecarts_types;

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
