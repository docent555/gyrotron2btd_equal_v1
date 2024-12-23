Lz = 5.243;
Lzi = 3.2;
Tend = 100;
Delta = 0.6;
Ia = 100;
R0 = 2.077;
Rb = 0.827;
g = 1.2;
ukv = 100;

Tend_ns = 1e-9 * Tend; % perevod v [ns]
gamma = 1.0 + ukv/511.0;
betta = sqrt(1.0d0 - 1.0/(gamma*gamma));
betta2 = 1.0d0 - 1.0/(gamma*gamma);
betta_z = betta/sqrt(g*g + 1.0d0);
betta_z2 = betta2/(g*g + 1.0d0);
betta_perp2 = betta2 - betta_z2;

c = 29979245800; % [cm/s]
e = 4.803e-10; % [ед. СГС]
m = 9.1093837015e-28; % [g]

% c = 299792458; % [m/s]
% e = 1.6021766208e-19; % [Кл]
% m = 9.1093837015e-31; % [кг]

% Вычислене nu (мода ТЕ28.12)
syms x;
dydx = matlabFunction(diff(besselj(28, x), x));
z=0:0.01:100;
% plot(z,dydx(z))
approx = 74;
nu = fzero(dydx, approx);
% nu = 73.952055635763557;

w_op = c * nu / R0;

I0 = 16 / (17045.81697831) * (Ia * betta_z * besselj(27,nu*Rb/R0)^2) / ...
        (gamma * betta_perp2^3 * (nu^2 - 28^2) * besselj(28, nu)^2);

ZetaEx = betta_perp2/2.0/betta_z*w_op*Lz/c;
ZetaExInter = betta_perp2/2.0/betta_z*w_op*Lzi/c;

TauEnd = betta_perp2^2*w_op*Tend_ns/8/betta_z2;



% Rr_file = load('d:\Alex\Documents\Work\novozhilova\22-09-23 Продольная структура и профиль рез-ра\RR2812.dat');
% fileID = fopen('RR2812.bin');
% if fileID < 0
%     fprintf('\nError of file open.\n');
%     pause;
% end
% Rr_file = fread(fileID,[100 2],'double');
% fclose(fileID);

% zax = z0:dz:zex;
% [ba, ca, da] = spline(length(Rr_file(:,1)), Rr_file(:,1), Rr_file(:,2));
%
% x = (zax(1:Nz3) - z0)/(zax(Nz3) - z0)*(Rr_file(end,1) - Rr_file(1,1)) + Rr_file(1,1);
% if (SPLINE == true)
%     for i=1:length(x)
%         Rr(i) = seval(x(i), length(Rr_file(:,1)), Rr_file(:,1), Rr_file(:,2), ba, ca, da);
%     end
%     Rr(:) = Rr(:)/10; % в сантиметрах
% else
%     for i=1:length(x)
%         Rr(i) = uval(x(i),Rr_file(:,1),Rr_file(:,2));
%     end
%     Rr(:) = Rr(:)/10; % в сантиметрах
% end

Rr_file = load('RR2812.dat');

Rr_file(:,2) = Rr_file(:,2)/10; % в сантиметрах

wc = c*nu*ones(length(Rr_file(:,2)),1)./Rr_file(:,2);    
Rr = 8*betta_z2/betta_perp2^2*(1 - wc/w_op);

dlt = [Rr_file(:,1) Rr];

save rr.dat dlt -ascii

% fileID = fopen("rr.dat");
% d = fscanf(fileID, '%f %f', [2 Inf]);
% fclose(fileID);
% d = d';
%plot(d(:,1), d(:,2))