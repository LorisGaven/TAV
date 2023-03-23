clear;
close all;

taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);


% Parametres :
R = 7;					% Rayon des disques
nb_points_affichage_disque = 30;
increment_angulaire = 2*pi/nb_points_affichage_disque;
theta = 0:increment_angulaire:2*pi;
rose = [253 108 158]/255;
q_max = 200;
nb_affichages = 1000;
pas_entre_affichages = floor(q_max/nb_affichages);
temps_pause = 0.0001;

beta = 1;
S = 140;
gamma = 5;
T_0 = 0.1;
lambda_0 = 100;
alpha = 0.99;
seuil = sqrt(2) * R;

% Lecture et affichage de l'image :
I = imread('colonie.png');
I = rgb2ycbcr(I);
I = double(I(:,:,1));
[nb_lignes,nb_colonnes] = size(I);
figure('Name',['Detection de flamants roses'],'Position',[0,0,L,0.58*H]);


lambda = lambda_0;
T = T_0;
c = [];
convergence = false;
q = 1;
while ~convergence

    % Naissances
    N = poissrnd(lambda);
    naissances = zeros(N,2);
    for i = 1:N
	    naissance_i = [nb_colonnes*rand nb_lignes*rand];
	    naissances(i,:) = naissance_i;
    end
    c = [c; naissances];

    % Tri des disques
    N_c = length(c);
    U_i = zeros(N_c, 1);
    for i = 1:N_c
        U_i(i) = 1 - (2 / (1 + exp(-gamma * ((calcul_I_moyen(I,c(i,:),R)/S) - 1))));
    end
    [U_i, idx] = sort(U_i, 'descend');
    c = c(idx, :);

    % Morts
    c_new = c;
    i = 1;
    while i <= N_c
        superpose = false;
        for j = 1:N_c
            if (j ~= i)
                superpose = superpose || norm(c_new(i,:) - c_new(j,:)) <= seuil;
            end
        end
        p = lambda / (lambda + exp((-U_i(i) - superpose) / T));
        if rand() < p
            c_new(i,:) = [];
            U_i(i,:) = [];
            N_c = N_c - 1;
        else
            i = i+1;
        end
    end

    % Convergence
    convergence = q == q_max;
    q = q + 1;
    T = alpha * T;
    lambda = alpha * lambda;
    c = c_new;

    % Affichage de la configuration initiale :
    subplot(1,2,1);
    imagesc(I);
    axis image;
    axis off;
    colormap gray;
    hold on;
    for i = 1:N_c
	    x_affich = c(i,1)+R*cos(theta);
	    y_affich = c(i,2)+R*sin(theta);
	    indices = find(x_affich>0 & x_affich<nb_colonnes & y_affich>0 & y_affich<nb_lignes);
	    plot(x_affich(indices),y_affich(indices),'Color',rose,'LineWidth',3);
    end
    pause(temps_pause);

end

fprintf('%d\n', N_c);
