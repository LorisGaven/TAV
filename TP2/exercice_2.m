donnees;
drawnow;

% Trace des donnees d'apprentissage (croix bleues) :
figure('Name','Estimation par les moindres carres','Position',[0.33*L,0,0.33*L,0.5*H]);
plot(D_app(1,:),D_app(2,:),'+b','MarkerSize',10,'LineWidth',2);
set(gca,'FontSize',20);
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
axis([-taille taille -taille taille]);
axis equal;
hold on;

% Tirages aleatoires de parametres pour l'ellipse :
nb_tirages = 10000;
parametres_test = zeros(nb_tirages,5);
parametres_test(:,1) = 2*taille/5*(rand(nb_tirages,1)+1);		% Demi-grand axe
parametres_test(:,2) = rand(nb_tirages,1);				% Excentricite
parametres_test(:,3) = (3*taille/5)*(2*rand(nb_tirages,1)-1);		% Abscisse du centre
parametres_test(:,4) = (3*taille/5)*(2*rand(nb_tirages,1)-1);		% Ordonnee du centre
parametres_test(:,5) = 2*pi*rand(nb_tirages,1);				% Angle du grand axe

% Estimation de l'ellipse par le maximum de vraisemblance :
parametres_MV = max_vraisemblance(D_app,parametres_test);

% Trace de l'ellipse estimee par le maximum de vraisemblance (trait rouge) :
[x_MV,y_MV] = points_ellipse(parametres_MV,theta_affichage);
plot([x_MV x_MV(1)],[y_MV y_MV(1)],'r-','LineWidth',3);

% Calcul et affichage du score :    
score_MV = calcul_score(parametres_VT,parametres_MV);
fprintf('Score de l''estimation par MV : %.3f\n',score_MV);

% Estimation en moindres carres :
X = moindres_carres(D_app);
parametres_MC = conversion(X);

% Trace de l'ellipse estimee en moindres carres (trait vert) :
[x_MC,y_MC] = points_ellipse(parametres_MC,theta_affichage);
plot([x_MC x_MC(1)],[y_MC y_MC(1)],'g-','LineWidth',3);

% Calcul et affichage du score :
score_MC = calcul_score(parametres_VT,parametres_MC);
fprintf('Score de l''estimation par MC : %.3f\n',score_MC);
legend(' Donnees d''apprentissage',' Ellipse estimee par MV',' Ellipse estimee par MC','Location','Best');

function mc = moindres_carres(D_app)
    n = size(D_app, 2);
    A = [D_app(1,:)'.^2 D_app(1,:)'.*D_app(2,:)' D_app(2,:)'.^2 D_app(1,:)' D_app(2,:)' ones(n,1)];
    A = [A; 1 0 1 0 0 0];
    B = [zeros(n,1); 1];
    mc = A \ B;
end

function mv = max_vraisemblance(D_app, param)
    nb_param = length(param);
    r_min = sum(calcul_r(D_app, param(1,:)).^2);
    mv = param(1,:);
    for i = 2:nb_param
        r = sum(calcul_r(D_app, param(i,:)).^2);
        if (r < r_min)
            r_min = r;
            mv = param(i,:);
        end
    end
end
