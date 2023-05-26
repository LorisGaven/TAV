clear;
close all;

% Paramètres à régler :
nb_iterations = 20;
epsilon = 0.01;
lambda = 100;

% Mise en place de la figure pour affichage :
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);
figure('Name','Decomposition structure + texture par le modele ROF','Position',[0,0,L,0.5*H]);

% Lecture et affichage de l'image u :
u = imread('Images/pilier.png');
u = double(u);
[nb_lignes,nb_colonnes,nb_canaux] = size(u);
subplot(1,3,1);
imagesc(uint8(u));
axis equal;
axis off;
title('Image','FontSize',30);

% Calcul des matrices Dx et Dy :
nb_pixels = nb_lignes*nb_colonnes;
e = ones(nb_pixels,1);
Dx = spdiags([-e e],[0 nb_lignes],nb_pixels,nb_pixels);
Dx(nb_pixels-nb_lignes+1:nb_pixels,:) = 0;
Dy = spdiags([-e e],[0 1],nb_pixels,nb_pixels);
Dy(nb_lignes:nb_lignes:nb_pixels,:) = 0;

% Seconds membres des systèmes à résoudre :
b = zeros(nb_lignes*nb_colonnes,nb_canaux);
for c = 1:nb_canaux
	u_c = u(:,:,c);
	b(:,c) = u_c(:);
end

% Affichage de la structure initiale :
u_barre = u;
subplot(1,3,2);
imagesc(uint8(u_barre));
axis equal;
axis off;
title('Structure (0%)','FontSize',30);

% Affichage de la texture initiale :
subplot(1,3,3);
imagesc(uint8(u-u_barre));
axis equal;
axis off;
title('Texture (0%)','FontSize',30);

drawnow;

for it = 1:nb_iterations

	for c = 1:nb_canaux

		% Mise à jour de la structure :
		u_barre(:,:,c) = calcul_structure_2(u_barre(:,:,c),b(:,c),Dx,Dy,lambda,epsilon);
	end
	
	% Affichage de la structure :
	subplot(1,3,2);
	imagesc(uint8(u_barre));
	axis equal;
	axis off;
	title(['Structure (' num2str(100*it/nb_iterations,'%.0f') '%)'],'FontSize',30);

	% Affichage de la texture :
	subplot(1,3,3);
	imagesc(uint8(u-u_barre));
	axis equal;
	axis off;
	title(['Texture (' num2str(100*it/nb_iterations,'%.0f') '%)'],'FontSize',30);

	drawnow;
end

function u_kp1 = calcul_structure_2(u_k,b,Dx,Dy,lambda,epsilon)
    [h, w] = size(u_k);
    u_k = u_k(:);
    n = length(u_k);
    grad_u_k = [Dx*u_k Dy*u_k];
    V = 1 ./ sqrt(grad_u_k(:,1).^2 + grad_u_k(:,2).^2 + epsilon);
    W = spdiags(V, 0, n, n);
    I = speye(n);
    A = I - lambda * (-Dx'*W*Dx - Dy'*W*Dy);
    u_kp1 = reshape(A\b, [h, w]);
end
