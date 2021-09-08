function [Psi0,error,X,Y,zX,zY,Pd,Ux,Uy]=MMP(V_in)%Multiple multipole method (MMP method) for Graphene
%clc
%clear all

hbar = 6.582119514e-13; %meV*s
c = 1e15; %nm/s
hbc = hbar*c; %meV*nm

phi =atan(0.5);
%parameters
R = 100; %nm

%Gate potential
h = 0.01/1;
% tk2R = 25:h:35;
% V2 = tk2R*hbc/R; %meV

%Incident energy
%E = 100*ones(1, length(V2)); %meV
%eta = 0.1; % eta = E/V2
%E = eta*V2;
%k1R = E*R/hbc; 

%k1R = 49.65;
%k1R = 69.9
%k1R = .1519
%k1R = 2.9537
%E = k1R*hbc/R;

% E = 1;
% k1R = E*R/hbc;
k1R = 0.1;
E = k1R*hbc/R;

V2 =V_in*hbc/R;

s2 = sign(E - V2);
k2R = abs(E - V2)*R/hbc;
tk2R = k1R;
%k2R = (0 + V2)*R/hbc;

%Geometry and Multipole setting
%set the order of multipoles
L = 1; % the maximum angular momentum
nL = -L:L;


%-----------------------Heart Shape---------------------------------
a = 0.3;
%set the # of multipoles ans assign their locations in polar coordinates
 Nm = 44; %inside the dot
 thetm = linspace(0, 2*pi, Nm);
 %Rm = 0.5; %Rm*R = rm
 Rm = 0.5*(1+a*cos(thetm));
 zRm = Rm.*exp(i*thetm);
%outside
 Nl = 48; %outside the dot
 thetl = linspace(0, 2*pi, Nl);
% %Rl = 2; %Rl*R = rl
 Rl = 2*(1+a*cos(thetl));
 zRl = Rl.*exp(i*thetl);
%discretization at the boundary: matching points #
 Nj = length(nL)*(Nm + Nl);
 theta = linspace(0, 2*pi, Nj);
%outward normal unit
 alpha = theta;
 %rho = 1; %rho*R = r
 rho = 1*(1+a*cos(theta));
 Z = rho.*exp(i*theta);
 X = real(Z); Y = imag(Z);
%------------------------Heart Shape-------------------------------

%Calculate the displacement matrix between the poles and matching points
[zm, zjm] = meshgrid(zRm, Z);
Djm = zjm - zm; 
Phjm = angle(Djm);
Rjm = abs(Djm);

[zl, zjl] = meshgrid(zRl, Z);
Djl = zjl - zl;
Phjl = angle(Djl);
Rjl = abs(Djl);


%incident wave
beta = 0*pi/4; %direction

%detection/measurement position
phi0 = pi/3;
r0 = .8;
z0 =  r0*exp(i*phi0); %for the stadium shape
%z0 = r0*exp(i*phi0);
z0l = z0 - zRl;
Ph0l = angle(z0l);
R0l = abs(z0l);

% %%the solving regions for the oval shape
Nx = 400;
Xxl =-1.5;
Xxu = 1.5;
%Xxl = .78;
%Xxu = 1.5;
Xx = linspace(Xxl, Xxu, Nx);
Yyl = -1.5;
Yyu = 1.5;
%Yyu = -.75;
Ny = floor(Nx*(Yyu - Yyl)/(Xxu - Xxl));
Yy = linspace(Yyl, Yyu, Ny);
[zX, zY] = meshgrid(Xx, Yy);
ZXY = zX + i*zY;
[Nr Nc] = size(ZXY);
Zxy = reshape(ZXY, 1, Nr*Nc);
%Ol = 0; Or = 0;

%Ind1 = find(real(Zxy) <= -a & abs(Zxy - Ol) <= 1);
%Oud1 = find(real(Zxy) <= -a & abs(Zxy - Ol) > 1);

%Ind2 = find(abs(real(Zxy)) < a & abs(imag(Zxy)) <= 1);
%Oud2 = find(abs(real(Zxy)) < a & abs(imag(Zxy)) > 1);

%Ind3 = find(real(Zxy) >= a & abs(Zxy - Or) <= 1);
%Oud3 = find(real(Zxy) >= a & abs(Zxy - Or) > 1);

%%the solving regions for rectangle shape
% rat = .5; Nx = 300;
% Xxl = -a - rat;
% Xxu = a + rat;
% Xx = linspace(Xxl, Xxu, Nx);
% Yyl = -1 - rat;
% Yyu = 1 + rat;
% Ny = floor(Nx*(Yyu - Yyl)/(Xxu - Xxl));
% Yy = linspace(Yyl, Yyu, Ny);
% [zX, zY] = meshgrid(Xx, Yy);
% ZXY = zX + i*zY;
% [Nr Nc] = size(ZXY)
% Zxy = reshape(ZXY, 1, Nr*Nc);
% %Ol = -a; Or = a;
% 
% Ind2 = find(abs(real(Zxy)) <= a & abs(imag(Zxy))<= 1);
% Oud2 = find(abs(real(Zxy)) > a | abs(imag(Zxy)) > 1);

%select set
Sin = ones(1, Nx*Ny); Sou = Sin;
Ind=find(inpolygon(real(Zxy),imag(Zxy),X,Y)==1);
Oud=find(inpolygon(real(Zxy),imag(Zxy),X,Y)==0);
tZin = Zxy; tZin(Oud) = 0; Sin(Oud) = 0;
tZou = Zxy; tZou(Ind) = Xxu; Sou(Ind) = 0;

tXou = real(tZou); tYou = imag(tZou);

[tzm, zOm] = meshgrid(zRm, tZou);
tDou = zOm - tzm; 
tPhou = angle(tDou);
tRou = abs(tDou);

[tzl, zIl] = meshgrid(zRl, tZin);
tDin = zIl - tzl;
tPhin = angle(tDin);
tRin = abs(tDin);

for n = 1:length(k2R)
    Ai = []; Bi = [];
    Ao = []; Bo = [];
    A0l = []; B0l=[]; C0l=[];
    Ain = []; Bin = []; Cin = [];
    Aou = []; Bou = []; Cou = [];
    InA = []; InB = []; InC = [];
    
    pInA = exp(i*k1R(n)*(cos(beta)*X' + sin(beta)*Y'));
    InB1 = cos(phi)*exp(-i*beta)*exp(i*k1R(n)*(cos(beta)*X' + sin(beta)*Y'));
    InB2 = sin(phi)*exp(i*beta)*exp(i*k1R(n)*(cos(beta)*X' + sin(beta)*Y'));
    pInB = cos(phi)*InB1.*exp(i*alpha') + sin(phi)*InB2.*exp(-i*alpha');
    pIn = [pInA; pInB];
    
    InA = exp(i*k1R(n)*(cos(beta)*tXou' + sin(beta)*tYou')).*Sou';
    InB = cos(phi)*exp(-i*beta)*exp(i*k1R(n)*(cos(beta)*tXou' + sin(beta)*tYou')).*Sou';
    InC = sin(phi)*exp(i*beta)*exp(i*k1R(n)*(cos(beta)*tXou' + sin(beta)*tYou')).*Sou';
    
    
    for l = 1:length(nL)
        
        HjmA = i*besselh(nL(l), 1, k1R(n)*Rjm).*exp(i*nL(l)*Phjm);
        HjmB1 = cos(phi)*besselh(nL(l)-1, 1, k1R(n)*Rjm).*exp(i*(nL(l)-1)*Phjm);
        HjmB2 = -sin(phi)*besselh(nL(l)+1, 1, k1R(n)*Rjm).*exp(i*(nL(l)+1)*Phjm);
        HjmB = cos(phi)*diag(exp(i*alpha))*HjmB1 + sin(phi)*diag(exp(-i*alpha))*HjmB2;
        
        %Ai(:,l) = sum(HjmA, 2);
        %Bi(:,l) = sum(HjmB, 2);
        Ai = [Ai HjmA];
        Bi = [Bi HjmB];
        
        HjlA = i*s2(n)*besselh(nL(l), 1, k2R(n)*Rjl).*exp(i*nL(l)*Phjl);
        HjlB1 = cos(phi)*besselh(nL(l)-1, 1, k2R(n)*Rjl).*exp(i*(nL(l)-1)*Phjl);
        HjlB2 = -sin(phi)*besselh(nL(l)+1, 1, k2R(n)*Rjl).*exp(i*(nL(l)+1)*Phjl);
        HjlB = cos(phi)*diag(exp(i*alpha))*HjlB1 + sin(phi)*diag(exp(-i*alpha))*HjlB2;
        
        %Ao(:,l) = sum(HjlA, 2);
        %Bo(:,l) = sum(HjlB, 2);    
        Ao = [Ao HjlA];
        Bo = [Bo HjlB];
        
        H0lA = i*s2(n)*besselh(nL(l), 1, k2R(n)*R0l).*exp(i*nL(l)*Ph0l);
        H0lB1 = cos(phi)*besselh(nL(l)-1, 1, k2R(n)*R0l).*exp(i*(nL(l)-1)*Ph0l);
        H0lB2 = -sin(phi)*besselh(nL(l)+1, 1, k2R(n)*R0l).*exp(i*(nL(l)+1)*Ph0l);
        %H0lB = H0lB1 + H0lB2;
        %A0l(l) = sum(H0lA);
        %B0l(l) = sum(H0lB);
        A0l = [A0l H0lA];
        B0l = [B0l H0lB1];
        C0l = [C0l H0lB2];
        
        HinA = i*s2(n)*besselh(nL(l), 1, k2R(n)*tRin).*exp(i*nL(l)*tPhin);
        HinB1 = cos(phi)*besselh(nL(l)-1, 1, k2R(n)*tRin).*exp(i*(nL(l)-1)*tPhin);
        HinB2 = -sin(phi)*besselh(nL(l)+1, 1, k2R(n)*tRin).*exp(i*(nL(l)+1)*tPhin);
        Ain = [Ain HinA];
        Bin = [Bin HinB1];
        Cin = [Cin HinB2];
        
        HouA = i*besselh(nL(l), 1, k1R(n)*tRou).*exp(i*nL(l)*tPhou);
        HouB1 = cos(phi)*besselh(nL(l)-1, 1, k1R(n)*tRou).*exp(i*(nL(l)-1)*tPhou);
        HouB2 = -sin(phi)*besselh(nL(l)+1, 1, k1R(n)*tRou).*exp(i*(nL(l)+1)*tPhou);
        Aou = [Aou HouA];
        Bou = [Bou HouB1];
        Cou = [Cou HouB2];
        
    end
    
    ABi = [Ai; Bi];  
    ABo = [Ao; Bo];
    [nr, nc] = size(ABi);
    M = [ABi, -ABo];
    
    C(:,n) = -pinv(M)*pIn;
    
    ee(n) = norm(M*C(:,n) + pIn)/norm(pIn);
    
    ABC0 = [A0l; B0l; C0l];
    psi0(:,n) = ABC0*C(nc+1:end, n);
    P0A = conj(psi0(1, n))*psi0(1, n); % psi2
    P0B = conj(psi0(2, n))*psi0(2, n); % psi1
    P0C = conj(psi0(3, n))*psi0(3, n); % psi3
    P0(n) = P0A + P0B + P0C;
    
    tPouA(:,n) = Aou*C(1:nc, n); PouA(:,n) = tPouA(:,n).*Sou' + InA;
    tPouB(:,n) = Bou*C(1:nc, n); PouB(:,n) = tPouB(:,n).*Sou' + InB;
    tPouC(:,n) = Cou*C(1:nc, n); PouC(:,n) = tPouC(:,n).*Sou' + InC;
    
    tPinA(:,n) = Ain*C(nc+1:end, n); PinA(:,n) = tPinA(:,n).*Sin';
    tPinB(:,n) = Bin*C(nc+1:end, n); PinB(:,n) = tPinB(:,n).*Sin';
    tPinC(:,n) = Cin*C(nc+1:end, n); PinC(:,n) = tPinC(:,n).*Sin';
    
    PsiA(:,n) = PinA(:,n) + PouA(:,n); % psi2
    PsiB(:,n) = PinB(:,n) + PouB(:,n); % psi1
    PsiC(:,n) = PinC(:,n) + PouC(:,n); % psi3
    
    PA(:,n) = conj(PsiA(:,n)).*PsiA(:,n);
    PB(:,n) = conj(PsiB(:,n)).*PsiB(:,n);
    PC(:,n) = conj(PsiC(:,n)).*PsiC(:,n);
    P(:,n) = PA(:,n) + PB(:,n) + PC(:,n);
end

Psi0 = P0 ;
error = ee;

%
Pn = P(:,n)'; 
Pd = reshape(Pn, Nr, Nc);
PdA = reshape(PA, Nr, Nc);
PdB = reshape(PB, Nr, Nc);
PdC = reshape(PC, Nr, Nc);

PsA = reshape(PsiA, Nr, Nc);
PsB = reshape(PsiB, Nr, Nc);
PsC = reshape(PsiC, Nr, Nc);

%current
Ux = 2*real(conj(PsA).*(cos(phi)*PsB + sin(phi)*PsC));
Uy = -2*imag(conj(PsA).*(cos(phi)*PsB - sin(phi)*PsC));
%Ux = Ux./sqrt(Ux.^2 + Uy.^2);
%Uy = Uy./sqrt(Ux.^2 + Uy.^2);



%Additional Plots
% indx = find(abs(Xx)<=1);
% xL = indx(1); xR = indx(end);
% zX(1, [xL xR])
% Pds = sum(Pd(:, xL:xR), 2);
% figure; 
% plot(Yy, log10(Pds))
% x = Yy;
% y = log10(Pds');
% lw=0.1/2; % width
% % plot
% lx=diff(x);
% ly=diff(y);
% lr=sqrt(lx.^2+ly.^2);
% lx=lx./lr;
% ly=ly./lr;
% Xe=x([1:end,end:-1:1])+lw*[ly(1), ly,-ly(end:-1:1), -ly(1)];
% Ye=y([1:end,end:-1:1])+lw*[-lx(1), -lx,lx(end:-1:1), lx(1)];
% figure; fill(Xe,Ye,'-r','EdgeColor','none','FaceAlpha',0.5)
% hold on; plot(Yy, log10(Pds))
% figure;
% plot(X, Y, 'k--.')
% hold on
% plot(real(zRm), imag(zRm), 'ro')
% plot(real(zRl), imag(zRl), 'b*')
% plot(real(z0), imag(z0), 'mp','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63],'MarkerSize',10)
% plot([real(Z); real(Zf)], [imag(Z); imag(Zf)], 'm')
% axis equal; 
end

