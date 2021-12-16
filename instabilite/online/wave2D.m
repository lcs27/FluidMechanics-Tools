function wave2D
%wave2D résoud la propagation d'ondes dans un milieu bipériodique
% pour différentes relations de dispersion:
%
% gravito-inertie: omega=sqrt(1+k2*coef)
%
% rossby         : omega=-kxx./(1+k2*coef)
% 
% interne        : omega=abs(kyy)./kk
%
% gravite        : omega=sqrt(kk.*tanh(kk*coef))
%
% gravite-courte : omega=sqrt(kk)
%
% gravite-longue : omega=kk
%
% où coef est un paramètre fixé
%
% Il calcule soit un problème aux conditions initiales,
% soit le sillage formé par un objet en déplacement,
% lui-même en translation ou oscillant


%% paramètres à régler par l'utilisateur
type_onde='gravite-longue';
type_probleme='CI'; % 'CI' ou 'objet'
type_mouvement='translation'; % 'translation' ou 'oscillation'
tend=40;%duree d'intégration
tplot=.5;%durée entre deux affichages à l'écran
coef=0.1; % coefficient entrant dans la relation de dispersion
damping=0e-4; % coefficient d'atténuation (dissipation)
U0=0.2; % vitesse de déplacement
theta=0; % angle du déplacement par rapport à l'axe des x (en degres)
nx=128;% nombre de points en x et en y (doit être une puissance de 2)

%% paramètres déduits + définition des variables
dx=2*pi/nx; % domaine [0 2*pi]x[0 2*pi]
x=(1:nx)*dx;% pas de grille
[xx,yy]=meshgrid(x,x);

kx=[0:nx/2 -nx/2+1:-1];% nombre d'onde selon x
[kxx,kyy]=meshgrid(kx,kx);
k2=kxx.^2+kyy.^2;
kk=sqrt(k2); % tableau contenant la norme du vecteur d'onde

dt=dx*.5; % pas de temps
nt=ceil(tend/dt); % nombre de pas de temps

r2=(xx-pi).^2+(yy-pi).^2;
sig=0.2; % largeur gaussienne
phi=exp(-r2/(2*sig^2)); % gaussienne de perturbation
hphi0=fft2(phi);%amplitude de Fourier complexe de la perturbation

% vecteur d'onde projeté selon le mouvement de l'objet
alpha=theta*pi/180;
kalpha=cos(alpha)*kxx+sin(alpha)*kyy;
if strcmp(type_probleme,'objet')
    hphi=hphi0*0;
else
    hphi=hphi0;
end

% on choisit la relation de dispersion
switch type_onde
    case 'gravito-inertie'
        omega=sqrt(1+k2*coef)*.2;
    case 'rossby'
        omega=-kxx./(1+k2*coef);
    case 'interne'
        omega=abs(kxx)./(kk+1e-6); % on ajoute 1e-6 pour éviter la division par 0
    case 'gravite'
        omega=sqrt(kk.*tanh(kk*coef))*.2;
    case 'gravite-courte'
        omega=sqrt(kk)*.5; 
    case 'gravite-longue'
        omega=kk*.5; % *.2 pour qu'elles n'aillent pas trop vite
end

% force à ce que la fenêtre soit à peu près carrée
clf
set(gcf,'position',[1 1 600 570])

maxi=max(abs(phi(:)));

%% boucle temporelle
for kt=0:nt
    t=kt*dt; %temps absolu
    
    % chaque amplitude est multipliée par son changement de phase
    % c'est le COEUR du calcul
    hphi=hphi.*exp(-1i*dt*omega);
    
    if strcmp(type_probleme,'objet')
        switch type_mouvement
            case 'translation'
                U=U0*(1-exp(-t/10)); % objet en translation
            case 'oscillation'
                U=U0*sin(t*2*pi/12); % object oscillant     
        end
        % l'objet se comporte comme une source pour l'amplitude complexe
        hphi=hphi-(dt*1i*U)*kalpha.*hphi0.*exp(-1i*t*U*kalpha);
%        hphi=hphi+dt*hphi0.*exp(-1i*t*U*kalpha);
    end
    
    if(damping~=0)
        % on ajoute du damping si besoin est
        hphi=hphi-hphi.*kk*damping;
    end
    
    if(floor(t/tplot)-floor( (t-dt)/tplot)==1) % n'affiche que tous les 'tplot'
        % reconstruit la solution dans l'espace physique=Fourier inverse
        phi=real(ifft2(hphi)); 
        pcolor(xx,yy,phi); shading flat
        axis xy;axis equal;
        %caxis([-1 1]*maxi*.2)
        xlabel('x')
        ylabel('y')
        title(sprintf('%s / time=%4.2f',type_onde,t))
        colorbar
        drawnow
    end
end