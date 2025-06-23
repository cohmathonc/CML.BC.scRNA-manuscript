
% This code is to generate full range 1-parameter saddle-node bifurcations for the scenario where
% individual cell types produce tristable response and the PsB also
% produces tristable response with the condition that sum of 'S' for all
% cell types is equal to the 'S' for the PsB.

% Notations
% C = BSR::ABL
% a11, a22, a33, a44 = ctPsBs
% U = transcriptome

clear
clc

Basin=[];
LocMinVar=[];

% parameters

kca1 = 200;
kua1 = 25;
gamma_a1 = 0.1653;
ca10 = 40;
ua10 = 400;
nca1 = 1;
nua1 = 6;

kca2 = 500;
kua2 = 22;
gamma_a2 = 0.1653;
ca20 = 100;
ua20 = 400;
nca2 = 1;
nua2 = 6;

kca3 = 550;
kua3 = 25;
gamma_a3 = 0.1653;
ca30 = 120;
ua30 = 400;
nca3 = 1;
nua3 = 6;

kca4 = 600;
kua4 = 28;
gamma_a4 = 0.1653;
ca40 = 300;
ua40 = 400;
nca4 = 1;
nua4 = 6;

ka1u = 40;
ka2u = 45;
ka3u = 48;
ka4u = 45;
a1u0 = 125;
a2u0 = 130;
a3u0 = 150;
a4u0 = 150;
na1u = 7;
na2u = 7;
na3u = 7;
na4u = 7;

kuu = 36.6837;
uu0 = 125;
nuu = 3;
gamma_u = 0.1653;



 
 w1 = 1;       % w1=1 for A1, w1=0 for No A1
 w2 = 1;       % w2=1 for A2, w2=0 for No A2
 w3 = 1;       % w3=1 for A3, w3=0 for No A3
 w4 = 1;       % w4=1 for A4, w4=0 for No A4



sig=[0:0.001:10];   % range of C



nr = length(sig);

U_SSS = NaN(nr,3);
U_USS = NaN(nr,3);
A_SSS = NaN(nr,3);
A_USS = NaN(nr,3);

cmap = jet(length(sig));

for ij=1:nr
s=sig(ij);

figure(1)
%subplot(5,5,1)



delb=01;

u1=[0:delb:3000];     % range for argument in potential function
        
        
        a11 = ((kca1.*(s./ca10).^nca1./(1+(s./ca10).^nca1)) ...
            +(kua1./(1+(u1./ua10).^nua1)))./gamma_a1;

        a21 = ((kca2.*(s./ca20).^nca2./(1+(s./ca20).^nca2)) ...
            +(kua2./(1+(u1./ua20).^nua2)))./gamma_a2;

        a31 = ((kca3.*(s./ca30).^nca3./(1+(s./ca30).^nca3)) ...
            +(kua3./(1+(u1./ua30).^nua3)))./gamma_a3;

        a41 = ((kca4.*(s./ca40).^nca4./(1+(s./ca40).^nca4)) ...
            +(kua4./(1+(u1./ua40).^nua4)))./gamma_a4;
        
        
        f =  (ka1u./(1+(a11./a1u0).^na1u)).*w1 ...
            +(ka2u./(1+(a21./a2u0).^na2u)).*w2...
            +(ka3u./(1+(a31./a3u0).^na3u)).*w3...
            +(ka4u./(1+(a41./a4u0).^na4u)).*w4 ...
            +(kuu.*(u1./uu0).^nuu./(1+(u1./uu0).^nuu))-gamma_u.*u1;
           
                
z=cumtrapz(-f.*delb);         % Potential is the negative of the function integrated
z1=z-min(z);

%a1 = a11;
 
TF_min=islocalmin(z1)';   % finding the local min
LocMinforVar = u1(TF_min==1);
LocMinforPot = z1(TF_min==1);
%LocMinforVar2 = a1(TF_min==1);


LocMin_Indx = find(TF_min==1);

TF_max=islocalmax(z1)';   % finding the local max
LocMaxforVar = u1(TF_max==1);
LocMaxforPot = z1(TF_max==1);
%LocMaxforVar2 = a1(TF_max==1)

LocMax_Indx = find(TF_max==1);

U_SSS(ij,1:length(LocMinforVar)) = LocMinforVar;
U_USS(ij,1:length(LocMaxforVar)) = LocMaxforVar;

% A_SSS(ij,1:length(LocMinforVar2)) = LocMinforVar2;
% A_USS(ij,1:length(LocMaxforVar2)) = LocMaxforVar2;



end

plot(sig,U_SSS,'black.')
hold on
plot(sig,U_USS,'r.')





