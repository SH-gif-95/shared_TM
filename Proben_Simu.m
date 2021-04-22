
clear all
clc
close all
global Nx Ny dx X Y k k0 KX KY KZ % für berechnungen nötige globals
global UD zD xD;
UD=[]; zD=[]; xD=[];%für data-dumping noetige globals

lambda= 0.830; k0=2*pi/lambda; k=k0; %mikrometer bzw mikrometer^{-1}
w0 = 28;                % strahltaille fuer function "LASER"
rayl= pi*w0.^2 /lambda; % rayleighlaenge fuer function "LASER"

Nx  = 200;
Ny  = 200;
dx  = lambda;
dy  = lambda;
dkx = 2*pi/Nx/dx;
dky = 2*pi/Ny/dy;
x   = dx* ([0:Nx-1]-Nx/2);
y   = dy* ([0:Ny-1]-Ny/2);
kx  = dkx* ([0:Nx-1]-(Nx)/2);
ky  = dky*  ([0:Ny-1]-(Ny)/2);
x=fftshift(x);
y=fftshift(y);
kx=fftshift(kx);
ky=fftshift(ky);

[X,Y] = meshgrid(x,y);
[KX,KY] = meshgrid(kx,ky);
KZ = sqrt(k.^2 -KX.^2 -KY.^2);

laser = LASER(X,Y,+rayl,rayl,k);

figure; imagesc(abs(fftshift(laser)))

n=1.5;
sigma=1;
sigma2=1;

Faktor1=2;
Faktor2=1;

N_tilts = 30*30;

%% Oberfläche rauer machen- nur die Amplitude der "großen" Störung wird verändert! 
%modellierung raue oberfläche als dünnes optisches element (phasenfaktor)

A_1 = rand(Nx,Ny);
A_1(1:Nx/2, 1:Ny/2) = A_1(1:Nx/2, 1:Ny/2)/1;
A_1(1:Nx/2, Ny/2:Ny) = A_1(1:Nx/2, Ny/2:Ny)/2;
A_1(Nx/2:Nx, 1:Ny/2) = A_1(Nx/2:Nx, 1:Ny/2)/3;
A_1(Nx/2:Nx, Ny/2:Ny) = A_1(Nx/2:Nx, Ny/2:Ny)/4;

p_1 = exp((-X.^2-Y.^2)./sigma^2);

surface_1 = ifft2( fft2(A_1) .* (fft2(p_1)));
surface_1 = Faktor1.*surface_1./max(surface_1(:));

A_2 = rand(Nx,Ny);
p_2 = exp((-X.^2-Y.^2)./sigma2^2);

surface_2=ifft2( fft2(A_2) .* (fft2(p_2)));
surface_2 = Faktor2.*surface_2./max(surface_2(:));

% surface_3 = surface_1 + surface_2;
surface_3 = surface_1;

% figure;imagesc(abs(surface3))
% figure;imagesc(angle(surface3))

probe_phi = 2*pi/lambda * ((n-1)*surface_3); %phasenverzögerung

figure; imagesc(probe_phi)

A=std(surface_3);
Rauheit=mean(A);

bsp = laser.*exp(1i*probe_phi);
grenz(n); %brechungsindexsprung

distance=1000;
ausgang = propagate(bsp,distance); %ausbreitung
figure; imagesc(abs(fftshift(ausgang)).^2);

%% scannender Laser

%% Grenzschicht
function dummy = grenz(n)
        global k k0 KX KY KZ
        k=n*k0;
        KZ = sqrt(k.^2 -KX.^2 -KY.^2);
end
%% free space propagation with 2 lateral axis: x,y
    function new = propagate(old,distance)
        global KZ
        if distance<0; %in the case of backwards propagation...
            KZ(imag(KZ)>0)=0; %...prevent that evanescent waves are amplified by the exp(1i*..)-term
        end
        %propagation equation with non-fftshifted axis
        something =  fft2(old) .* exp(1i*KZ.*distance);
        new =  ifft2(something);
    end
%% Drehung des Schirms
    function neu = drehung(alt,winkelgrad)
        global k KX KY KZ Nx Ny
        %disp(sprintf('Drehung des Schirms um %d Grad',winkelgrad))
        
        alpha=winkelgrad*pi/180;

        K_xyz=[reshape(KX,Nx*Ny,1) , reshape(KY,Nx*Ny,1)  , reshape(KZ,Nx*Ny,1) ];
        %rotX =[1,0,0 ; 0, cos(alpha),-sin(alpha);  0,sin(alpha),cos(alpha)];
        rotY =[cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];
        %rotZ =[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
        K_abc  = K_xyz*rotY;

        K_abcCOPY=K_abc;
        K_abc(imag(K_xyz(:,3)) ~= 0,:) = 12323;
        K_a=reshape(K_abc(:,1),Ny,Nx);
        K_b=reshape(K_abc(:,2),Ny,Nx);
        K_c=reshape(K_abc(:,3),Ny,Nx);

        K_aC=reshape(K_abcCOPY(:,1),Ny,Nx);
        K_bC=reshape(K_abcCOPY(:,2),Ny,Nx);
        K_cC=reshape(K_abcCOPY(:,3),Ny,Nx);

        A=(fft2( (alt) ));
        A(K_aC.^2 +K_bC.^2 >= k.^2) = 0;

        %interp2 kann nur mit monotonen arrays durchgefuehrtwerden->shiften
        AI = interp2(fftshift(KX),fftshift(KY),fftshift(A),fftshift(K_a),fftshift(K_b),'linear');
        AI = fftshift(AI); 
        AI(KX.^2 +KY.^2 > k.^2) = 0.0;
        
        AI=AI.* K_c./KZ;
        AI(isnan(AI))=0;
        AI(abs(AI)==Inf)=0;
       
        neu=ifft2(AI);
    end
%% GAUSSSTRAHL - Funktionen
    function [invRz] = invR(z,zr)
        invRz = z./(z.^2 + zr^2);
    end
    function [Wz] = W(z,zr,k)
        w0 = sqrt(zr*2/k);
        Wz = w0*sqrt( 1 + (z./zr).^2 );
    end
    function [A]=LASER(x,y,z,zr,k)
        disp(sprintf('erzeuge Gausstrahl...'))
        w0 = sqrt(zr*2/k);
        A = 1 * w0./W(z,zr,k) .* exp( -(x.^2+y.^2) ./(W(z,zr,k).^2) - j*( k.*(x.^2+y.^2).*invR(z,zr)/2 - atan2(z,zr) )); 
    
    end
    function [A] = easygauss(x,y,zr,k)
        w0 = sqrt(zr*2/k);
        A = exp(-(x.^2+y.^2) ./ w0^2);
    end
    function [rect] = rect(x)
    rect=0.*x+1;
    rect(find(x>0.5))=0;
    rect(find(x<-0.5))=0;    
    
    end
%% zusaetzliche Datenspeicherung (Zentrum, Z-Position, Querschnitt....)
    function DUMP(feld,PosZ)
        global Ny Nx dx
        global UD xD zD
       
        [dummy,indx] = max(feld(Nx/2+1,:));
        
        xD=[xD, indx*dx];
        zD=[zD, PosZ];
        UD = [UD;feld(floor(Nx/2),:)];
    end

