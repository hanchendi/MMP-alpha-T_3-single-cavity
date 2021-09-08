%Multiple multipole method (MMP method) for Graphene
clc
clear all

hbar = 6.582119514e-13; %meV*s
c = 1e15; %nm/s
hbc = hbar*c; %meV*nm

phi = atan(0.5);
%parameters
R = 100; %nm

%Gate potential
h = 0.002;
% tk2R = 25:h:35;
% V2 = tk2R*hbc/R; %meV

V2R = 10:h:12;
V2 = V2R*hbc/R;
%ks = (.5 - .02)*V2R;
%kf = (.5 + .01)*V2R;
%k1R = ks:h:kf;
%k1R = 27:h:29;
k1R =0.1*ones(1,length(V2R)) ;
E = k1R*hbc/R;

s2 = sign(E - V2);
k2R = abs(E - V2)*R/hbc;
tk2R = k1R;

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
%------------------------Eplltic Shape-------------------------------

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
% pInA = exp(-i*beta/2)*exp(i*k1R*(cos(beta)*X' + sin(beta)*Y'));
% pInB = exp(i*beta/2)*exp(i*k1R*(cos(beta)*X' + sin(beta)*Y'));
% pIn = [pInA; pInB];

%detection/measurement position
phi0 = pi/3;
r0 = .8;
z0 =  r0*exp(i*phi0); %for the stadium shape
%z0 = r0*exp(i*phi0);
z0l = z0 - zRl;
Ph0l = angle(z0l);
R0l = abs(z0l);


for n = 1:length(V2R)
    Ai = []; Bi = [];
    Ao = []; Bo = [];
    A0l = []; B0l=[];
   
    
    pInA = exp(i*k1R(n)*(cos(beta)*X' + sin(beta)*Y'));
    InB1 = cos(phi)*exp(-i*beta)*exp(i*k1R(n)*(cos(beta)*X' + sin(beta)*Y'));
    InB2 = sin(phi)*exp(i*beta)*exp(i*k1R(n)*(cos(beta)*X' + sin(beta)*Y'));
    pInB = cos(phi)*InB1.*exp(i*alpha') + sin(phi)*InB2.*exp(-i*alpha');
    pIn = [pInA; pInB];
    
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
        H0lB = H0lB1 + H0lB2;
        %A0l(l) = sum(H0lA);
        %B0l(l) = sum(H0lB);
        A0l = [A0l H0lA];
        B0l = [B0l H0lB];
    end
    
    ABi = [Ai; Bi];  
    ABo = [Ao; Bo];
    [nr, nc] = size(ABi);
    M = [ABi, -ABo];
    
    C = -pinv(M)*pIn;
    
    ee(n) = norm(M*C + pIn)/norm(pIn);
    
    AB0 = [A0l; B0l];
    psi0(:) = AB0*C(nc+1:end);
    P0A = conj(psi0(1))*psi0(1);
    P0B = conj(psi0(2))*psi0(2);
    P0(n) = P0A + P0B;
    disp(n/length(V2R))
end

save data_heart_a_03