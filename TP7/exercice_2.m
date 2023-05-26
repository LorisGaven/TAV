clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

% Lecture et affichage de l'image source s :
figure('Name','Photomontage naif','Position',[0.1*L,0.1*H,0.9*L,0.7*H]);
s = imread('Images/rose.jpg');
[nb_lignes_s,nb_colonnes_s,nb_canaux] = size(s);
subplot(1,2,1);
imagesc(s);
axis image off;
title('Image source','FontSize',20);
hold on;

% Selection et affichage d'un polygone p dans s :
disp('Selectionnez un polygone (double-clic pour valider)');
[p,x_p,y_p] = roipoly(s);
for k = 1:length(x_p)-1
	line([x_p(k) x_p(k+1)],[y_p(k) y_p(k+1)],'Color','r','LineWidth',2);
end

% Bornes du rectangle englobant de p :
i_p = min(max(round(y_p),1),nb_lignes_s);
j_p = min(max(round(x_p),1),nb_colonnes_s);
i_p_min = min(i_p(:));
i_p_max = max(i_p(:));
j_p_min = min(j_p(:));
j_p_max = max(j_p(:));

% Lecture et affichage de l'image cible c :
c = imread('Images/rose.jpg');
c = rgb2lab(c);
[nb_lignes_c,nb_colonnes_c,nb_canaux] = size(c);
subplot(1,2,2);
imagesc(c);
axis image off;
title('Image cible','FontSize',20);
hold on;

% Selection et affichage d'un rectangle r dans c :
disp('Cliquez les deux extremites de la zone cible');
[x_r,y_r] = ginput(2);
i_r = min(max(round(y_r),1),nb_lignes_c);
j_r = min(max(round(x_r),1),nb_colonnes_c);
j_r_min = min(j_r(:));
j_r_max = max(j_r(:));
i_r_min = min(i_r(:));
i_r_max = max(i_r(:));
line([j_r_min j_r_max],[i_r_min,i_r_min],'Color','r','LineWidth',2);
line([j_r_min j_r_max],[i_r_max,i_r_max],'Color','r','LineWidth',2);
line([j_r_min j_r_min],[i_r_min,i_r_max],'Color','r','LineWidth',2);
line([j_r_max j_r_max],[i_r_min,i_r_max],'Color','r','LineWidth',2);

% Sous-matrice de c correspondant au rectangle r :
r = c(i_r_min:i_r_max,j_r_min:j_r_max,:);

% Seules les sous-matrices a l'interieur du rectangle englobant de p sont conservees :
s = s(i_p_min:i_p_max,j_p_min:j_p_max,:);
p = p(i_p_min:i_p_max,j_p_min:j_p_max);

% Redimensionnement de s et p aux dimensions de r :
[nb_lignes_r,nb_colonnes_r,nb_canaux] = size(r);
s = imresize(s,[nb_lignes_r,nb_colonnes_r]);
p = imresize(p,[nb_lignes_r,nb_colonnes_r]);

% Operateur gradient :
nb_pixels = nb_lignes_r * nb_colonnes_r;
e = ones(nb_pixels,1);
Dx = spdiags([-e e],[0 nb_lignes_r],nb_pixels,nb_pixels);
Dx(end-nb_lignes_r+1:end,:) = 0;
Dy = spdiags([-e e],[0 1],nb_pixels,nb_pixels);
Dy(nb_lignes_r:nb_lignes_r:end,:) = 0;

% Calcul et affichage de l'image resultat u :
u = c;
interieur = find(p>0);
u(i_r_min:i_r_max,j_r_min:j_r_max,:) = collage(r,s,interieur, Dx, Dy);
hold off;
imagesc(u);
axis image off;
title('Resultat du photomontage','FontSize',20);

function u = collage(r,s,interieur, Dx, Dy)
    [h,w,c] = size(r);
    r = double(r);
    s = double(s);
    u = zeros(h,w,c);
    bord_r = zeros(h,w);
    bord_r(1,:) = 1;
    bord_r(h,:) = 1;
    bord_r(:,1) = 1;
    bord_r(:,w) = 1;
    indices_bord_r = find(bord_r == 1);
    n_bord_r = length(indices_bord_r);
    n_r = h*w;

    A = - Dx'*Dx - Dy'*Dy;
    A(indices_bord_r,:) = sparse(1:n_bord_r,indices_bord_r,ones(n_bord_r,1),n_bord_r,n_r);

    for k = 1:size(r,3)
        rk = r(:,:,k);
        sk = s(:,:,k);
        dr_x = Dx*rk(:);
        dr_y = Dy*rk(:);
        grad_r = [dr_x dr_y];
        ds_x = Dx*sk(:);
        ds_y = Dy*sk(:);
        grad_s = [ds_x ds_y];
    
        g = grad_r;
        g(interieur,:) = grad_s(interieur,:);
    
        b = Dx*g(:,1) + Dy*g(:,2);
        b(indices_bord_r) = rk(indices_bord_r);

        u(:,:,k) = reshape(A\b, [h, w]);
    end
end
