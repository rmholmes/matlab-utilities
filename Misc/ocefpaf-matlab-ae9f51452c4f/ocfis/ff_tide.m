function [Amp, Pha, Tid] = ff_tide(data)
% ff_tide.m -> Harmonic analysis using least-squares fit
%
% use:  [Amp, Pha, Tid] = ff_tide(data)
% input:
%          data -> Sea Level High
%
% example:
%          load bl_2001
%          [Amp, Pha, Tid] = ff_tide(WL) % Water Level
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     29-May-2009
% modified: 16-Nov-2010
%
% obs: Decribed in Schuremann, 1941 manual of harmonic
%      analysis and prediction of tides. Washington, D.C.,U.S.
%      Coast & Geodetic Surv., S.P. n. 98, 317p.
%
%  Must be hourly series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Tidal analysis by least squares                    %
%                 Method from Schuremann (Schuremann, 1941)               %
%                                                                         %
%    SCHUREMANN, P. 1941.                                                 %
%    Manual of harmonic analysis and prediction of tides.                 %
%    Washington, D.C.,U.S. Coast & Geodetic Surv., S.P. n. 98, 317p.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Componentes de Longo Período
Sa  = 0.000717; Ssa= 0.001434;
Mm  = 0.009501; Mf = 0.019164;
Msf = 0.017731;

% Componentes de Maré Diurnas
Q1  = 0.233851; O1 = 0.243352;
M1  = 0.252934; P1 = 0.261083;
K1  = 0.262516; J1 = 0.272017;
OO1 = 0.281680;

% Componentes de Maré Semi-Diurnas
MNS2 = 0.478636; N22 = 0.486866;
MU2  = 0.488137; N2  = 0.496367;
NU2  = 0.497638; M2  = 0.505868;
L2   = 0.515369; T2  = 0.522882;
S2   = 0.523599; K2  = 0.525032;

% Componentes de Maré Ter-diurnas
MO3 = 0.749220; M3  = 0.758802;
MK3 = 0.768384;

% Componentes de Maré Quarto-diurnas
MN4 = 1.002235; M4  = 1.011736;
SN4 = 1.019966; MS4 = 1.029467;

%% Frequência Angular(rad/h)
% Componentes de Longo Período
w(1)  = Sa;  w(3)  = Ssa;
w(5)  = Mm;  w(7)  = Mf;
w(9)  = Msf;
% Componentes de Maré Diurnas
w(11) = Q1; w(13) = O1;
w(15) = M1; w(17) = P1;
w(19) = K1; w(21) = J1;
w(23) = OO1;
% Componentes de Maré Semi-Diurnas
w(25) = MNS2; w(27) = N22;
w(29) = MU2;  w(31) = N2;
w(33) = NU2;  w(35) = M2;
w(37) = L2;   w(39) = T2;
w(41) = S2;   w(43) = K2;
% Componentes de Maré Ter-diurnas
w(45) = MO3; w(47) = M3;
w(49) = MK3;
% Componentes de Maré Quarto-diurnas
w(51) = MN4; w(53) = M4;
w(55) = SN4; w(57) = MS4;

nhrs = max( size(data) );
data = data - mean(data);

tam = length(w)+1;

MT = zeros(tam,tam);
XT = zeros(1,tam);

for i =1:2:tam-1
  for t=1:nhrs
    MT(i,1) =   MT(i,1)   + cos(Sa  *t)*cos(w(i)*t);
    MT(i,2) =   MT(i,2)   + sin(Sa  *t)*cos(w(i)*t);
    MT(i,3) =   MT(i,3)   + cos(Ssa *t)*cos(w(i)*t);
    MT(i,4) =   MT(i,4)   + sin(Ssa *t)*cos(w(i)*t);
    MT(i,5) =   MT(i,5)   + cos(Mm  *t)*cos(w(i)*t);
    MT(i,6) =   MT(i,6)   + sin(Mm  *t)*cos(w(i)*t);
    MT(i,7) =   MT(i,7)   + cos(Mf  *t)*cos(w(i)*t);
    MT(i,8) =   MT(i,8)   + sin(Mf  *t)*cos(w(i)*t);
    MT(i,9) =   MT(i,9)   + cos(Msf *t)*cos(w(i)*t);
    MT(i,10)=   MT(i,10)  + sin(Msf *t)*cos(w(i)*t);
    MT(i,11)=   MT(i,11)  + cos(Q1  *t)*cos(w(i)*t);
    MT(i,12)=   MT(i,12)  + sin(Q1  *t)*cos(w(i)*t);
    MT(i,13)=   MT(i,13)  + cos(O1  *t)*cos(w(i)*t);
    MT(i,14)=   MT(i,14)  + sin(O1  *t)*cos(w(i)*t);
    MT(i,15)=   MT(i,15)  + cos(M1  *t)*cos(w(i)*t);
    MT(i,16)=   MT(i,16)  + sin(M1  *t)*cos(w(i)*t);
    MT(i,17)=   MT(i,17)  + cos(P1  *t)*cos(w(i)*t);
    MT(i,18)=   MT(i,18)  + sin(P1  *t)*cos(w(i)*t);
    MT(i,19)=   MT(i,19)  + cos(K1  *t)*cos(w(i)*t);
    MT(i,20)=   MT(i,20)  + sin(K1  *t)*cos(w(i)*t);
    MT(i,21)=   MT(i,21)  + cos(J1  *t)*cos(w(i)*t);
    MT(i,22)=   MT(i,22)  + sin(J1  *t)*cos(w(i)*t);
    MT(i,23)=   MT(i,23)  + cos(OO1 *t)*cos(w(i)*t);
    MT(i,24)=   MT(i,24)  + sin(OO1 *t)*cos(w(i)*t);
    MT(i,25)=   MT(i,25)  + cos(MNS2*t)*cos(w(i)*t);
    MT(i,26)=   MT(i,26)  + sin(MNS2*t)*cos(w(i)*t);
    MT(i,27)=   MT(i,27)  + cos(N22 *t)*cos(w(i)*t);
    MT(i,28)=   MT(i,28)  + sin(N22 *t)*cos(w(i)*t);
    MT(i,29)=   MT(i,29)  + cos(MU2 *t)*cos(w(i)*t);
    MT(i,30)=   MT(i,30)  + sin(MU2 *t)*cos(w(i)*t);
    MT(i,31)=   MT(i,31)  + cos(N2  *t)*cos(w(i)*t);
    MT(i,32)=   MT(i,32)  + sin(N2  *t)*cos(w(i)*t);
    MT(i,33)=   MT(i,33)  + cos(NU2 *t)*cos(w(i)*t);
    MT(i,34)=   MT(i,34)  + sin(NU2 *t)*cos(w(i)*t);
    MT(i,35)=   MT(i,35)  + cos(M2  *t)*cos(w(i)*t);
    MT(i,36)=   MT(i,36)  + sin(M2  *t)*cos(w(i)*t);
    MT(i,37)=   MT(i,37)  + cos(L2  *t)*cos(w(i)*t);
    MT(i,38)=   MT(i,38)  + sin(L2  *t)*cos(w(i)*t);
    MT(i,39)=   MT(i,39)  + cos(T2  *t)*cos(w(i)*t);
    MT(i,40)=   MT(i,40)  + sin(T2  *t)*cos(w(i)*t);
    MT(i,41)=   MT(i,41)  + cos(S2  *t)*cos(w(i)*t);
    MT(i,42)=   MT(i,42)  + sin(S2  *t)*cos(w(i)*t);
    MT(i,43)=   MT(i,43)  + cos(K2  *t)*cos(w(i)*t);
    MT(i,44)=   MT(i,44)  + sin(K2  *t)*cos(w(i)*t);
    MT(i,45)=   MT(i,45)  + cos(MO3 *t)*cos(w(i)*t);
    MT(i,46)=   MT(i,46)  + sin(MO3 *t)*cos(w(i)*t);
    MT(i,47)=   MT(i,47)  + cos(M3  *t)*cos(w(i)*t);
    MT(i,48)=   MT(i,48)  + sin(M3  *t)*cos(w(i)*t);
    MT(i,49)=   MT(i,49)  + cos(MK3 *t)*cos(w(i)*t);
    MT(i,50)=   MT(i,50)  + sin(MK3 *t)*cos(w(i)*t);
    MT(i,51)=   MT(i,51)  + cos(MN4 *t)*cos(w(i)*t);
    MT(i,52)=   MT(i,52)  + sin(MN4 *t)*cos(w(i)*t);
    MT(i,53)=   MT(i,53)  + cos(M4  *t)*cos(w(i)*t);
    MT(i,54)=   MT(i,54)  + sin(M4  *t)*cos(w(i)*t);
    MT(i,55)=   MT(i,55)  + cos(SN4 *t)*cos(w(i)*t);
    MT(i,56)=   MT(i,56)  + sin(SN4 *t)*cos(w(i)*t);
    MT(i,57)=   MT(i,57)  + cos(MS4 *t)*cos(w(i)*t);
    MT(i,58)=   MT(i,58)  + sin(MS4 *t)*cos(w(i)*t);

    MT(i+1,1) =   MT(i+1,1)   + cos(Sa  *t)*sin(w(i)*t);
    MT(i+1,2) =   MT(i+1,2)   + sin(Sa  *t)*sin(w(i)*t);
    MT(i+1,3) =   MT(i+1,3)   + cos(Ssa *t)*sin(w(i)*t);
    MT(i+1,4) =   MT(i+1,4)   + sin(Ssa *t)*sin(w(i)*t);
    MT(i+1,5) =   MT(i+1,5)   + cos(Mm  *t)*sin(w(i)*t);
    MT(i+1,6) =   MT(i+1,6)   + sin(Mm  *t)*sin(w(i)*t);
    MT(i+1,7) =   MT(i+1,7)   + cos(Mf  *t)*sin(w(i)*t);
    MT(i+1,8) =   MT(i+1,8)   + sin(Mf  *t)*sin(w(i)*t);
    MT(i+1,9) =   MT(i+1,9)   + cos(Msf *t)*sin(w(i)*t);
    MT(i+1,10)=   MT(i+1,10)  + sin(Msf *t)*sin(w(i)*t);
    MT(i+1,11)=   MT(i+1,11)  + cos(Q1  *t)*sin(w(i)*t);
    MT(i+1,12)=   MT(i+1,12)  + sin(Q1  *t)*sin(w(i)*t);
    MT(i+1,13)=   MT(i+1,13)  + cos(O1  *t)*sin(w(i)*t);
    MT(i+1,14)=   MT(i+1,14)  + sin(O1  *t)*sin(w(i)*t);
    MT(i+1,15)=   MT(i+1,15)  + cos(M1  *t)*sin(w(i)*t);
    MT(i+1,16)=   MT(i+1,16)  + sin(M1  *t)*sin(w(i)*t);
    MT(i+1,17)=   MT(i+1,17)  + cos(P1  *t)*sin(w(i)*t);
    MT(i+1,18)=   MT(i+1,18)  + sin(P1  *t)*sin(w(i)*t);
    MT(i+1,19)=   MT(i+1,19)  + cos(K1  *t)*sin(w(i)*t);
    MT(i+1,20)=   MT(i+1,20)  + sin(K1  *t)*sin(w(i)*t);
    MT(i+1,21)=   MT(i+1,21)  + cos(J1  *t)*sin(w(i)*t);
    MT(i+1,22)=   MT(i+1,22)  + sin(J1  *t)*sin(w(i)*t);
    MT(i+1,23)=   MT(i+1,23)  + cos(OO1 *t)*sin(w(i)*t);
    MT(i+1,24)=   MT(i+1,24)  + sin(OO1 *t)*sin(w(i)*t);
    MT(i+1,25)=   MT(i+1,25)  + cos(MNS2*t)*sin(w(i)*t);
    MT(i+1,26)=   MT(i+1,26)  + sin(MNS2*t)*sin(w(i)*t);
    MT(i+1,27)=   MT(i+1,27)  + cos(N22 *t)*sin(w(i)*t);
    MT(i+1,28)=   MT(i+1,28)  + sin(N22 *t)*sin(w(i)*t);
    MT(i+1,29)=   MT(i+1,29)  + cos(MU2 *t)*sin(w(i)*t);
    MT(i+1,30)=   MT(i+1,30)  + sin(MU2 *t)*sin(w(i)*t);
    MT(i+1,31)=   MT(i+1,31)  + cos(N2  *t)*sin(w(i)*t);
    MT(i+1,32)=   MT(i+1,32)  + sin(N2  *t)*sin(w(i)*t);
    MT(i+1,33)=   MT(i+1,33)  + cos(NU2 *t)*sin(w(i)*t);
    MT(i+1,34)=   MT(i+1,34)  + sin(NU2 *t)*sin(w(i)*t);
    MT(i+1,35)=   MT(i+1,35)  + cos(M2  *t)*sin(w(i)*t);
    MT(i+1,36)=   MT(i+1,36)  + sin(M2  *t)*sin(w(i)*t);
    MT(i+1,37)=   MT(i+1,37)  + cos(L2  *t)*sin(w(i)*t);
    MT(i+1,38)=   MT(i+1,38)  + sin(L2  *t)*sin(w(i)*t);
    MT(i+1,39)=   MT(i+1,39)  + cos(T2  *t)*sin(w(i)*t);
    MT(i+1,40)=   MT(i+1,40)  + sin(T2  *t)*sin(w(i)*t);
    MT(i+1,41)=   MT(i+1,41)  + cos(S2  *t)*sin(w(i)*t);
    MT(i+1,42)=   MT(i+1,42)  + sin(S2  *t)*sin(w(i)*t);
    MT(i+1,43)=   MT(i+1,43)  + cos(K2  *t)*sin(w(i)*t);
    MT(i+1,44)=   MT(i+1,44)  + sin(K2  *t)*sin(w(i)*t);
    MT(i+1,45)=   MT(i+1,45)  + cos(MO3 *t)*sin(w(i)*t);
    MT(i+1,46)=   MT(i+1,46)  + sin(MO3 *t)*sin(w(i)*t);
    MT(i+1,47)=   MT(i+1,47)  + cos(M3  *t)*sin(w(i)*t);
    MT(i+1,48)=   MT(i+1,48)  + sin(M3  *t)*sin(w(i)*t);
    MT(i+1,49)=   MT(i+1,49)  + cos(MK3 *t)*sin(w(i)*t);
    MT(i+1,50)=   MT(i+1,50)  + sin(MK3 *t)*sin(w(i)*t);
    MT(i+1,51)=   MT(i+1,51)  + cos(MN4 *t)*sin(w(i)*t);
    MT(i+1,52)=   MT(i+1,52)  + sin(MN4 *t)*sin(w(i)*t);
    MT(i+1,53)=   MT(i+1,53)  + cos(M4  *t)*sin(w(i)*t);
    MT(i+1,54)=   MT(i+1,54)  + sin(M4  *t)*sin(w(i)*t);
    MT(i+1,55)=   MT(i+1,55)  + cos(SN4 *t)*sin(w(i)*t);
    MT(i+1,56)=   MT(i+1,56)  + sin(SN4 *t)*sin(w(i)*t);
    MT(i+1,57)=   MT(i+1,57)  + cos(MS4 *t)*sin(w(i)*t);
    MT(i+1,58)=   MT(i+1,58)  + sin(MS4 *t)*sin(w(i)*t);
  end
end

for t=1:nhrs
  XT(1) =XT(1)  + data(t)*cos(Sa  *t);
  XT(2) =XT(2)  + data(t)*sin(Sa  *t);
  XT(3) =XT(3)  + data(t)*cos(Ssa *t);
  XT(4) =XT(4)  + data(t)*sin(Ssa *t);
  XT(5) =XT(5)  + data(t)*cos(Mm  *t);
  XT(6) =XT(6)  + data(t)*sin(Mm  *t);
  XT(7) =XT(7)  + data(t)*cos(Mf  *t);
  XT(8) =XT(8)  + data(t)*sin(Mf  *t);
  XT(9) =XT(9)  + data(t)*cos(Msf *t);
  XT(10)=XT(10) + data(t)*sin(Msf *t);
  XT(11)=XT(11) + data(t)*cos(Q1  *t);
  XT(12)=XT(12) + data(t)*sin(Q1  *t);
  XT(13)=XT(13) + data(t)*cos(O1  *t);
  XT(14)=XT(14) + data(t)*sin(O1  *t);
  XT(15)=XT(15) + data(t)*cos(M1  *t);
  XT(16)=XT(16) + data(t)*sin(M1  *t);
  XT(17)=XT(17) + data(t)*cos(P1  *t);
  XT(18)=XT(18) + data(t)*sin(P1  *t);
  XT(19)=XT(19) + data(t)*cos(K1  *t);
  XT(20)=XT(20) + data(t)*sin(K1  *t);
  XT(21)=XT(21) + data(t)*cos(J1  *t);
  XT(22)=XT(22) + data(t)*sin(J1  *t);
  XT(23)=XT(23) + data(t)*cos(OO1 *t);
  XT(24)=XT(24) + data(t)*sin(OO1 *t);
  XT(25)=XT(25) + data(t)*cos(MNS2*t);
  XT(26)=XT(26) + data(t)*sin(MNS2*t);
  XT(27)=XT(27) + data(t)*cos(N22 *t);
  XT(28)=XT(28) + data(t)*sin(N22 *t);
  XT(29)=XT(29) + data(t)*cos(MU2 *t);
  XT(30)=XT(30) + data(t)*sin(MU2 *t);
  XT(31)=XT(31) + data(t)*cos(N2  *t);
  XT(32)=XT(32) + data(t)*sin(N2  *t);
  XT(33)=XT(33) + data(t)*cos(NU2 *t);
  XT(34)=XT(34) + data(t)*sin(NU2 *t);
  XT(35)=XT(35) + data(t)*cos(M2  *t);
  XT(36)=XT(36) + data(t)*sin(M2  *t);
  XT(37)=XT(37) + data(t)*cos(L2  *t);
  XT(38)=XT(38) + data(t)*sin(L2  *t);
  XT(39)=XT(39) + data(t)*cos(T2  *t);
  XT(40)=XT(40) + data(t)*sin(T2  *t);
  XT(41)=XT(41) + data(t)*cos(S2  *t);
  XT(42)=XT(42) + data(t)*sin(S2  *t);
  XT(43)=XT(43) + data(t)*cos(K2  *t);
  XT(44)=XT(44) + data(t)*sin(K2  *t);
  XT(45)=XT(45) + data(t)*cos(MO3 *t);
  XT(46)=XT(46) + data(t)*sin(MO3 *t);
  XT(47)=XT(47) + data(t)*cos(M3  *t);
  XT(48)=XT(48) + data(t)*sin(M3  *t);
  XT(49)=XT(49) + data(t)*cos(MK3 *t);
  XT(50)=XT(50) + data(t)*sin(MK3 *t);
  XT(51)=XT(51) + data(t)*cos(MN4 *t);
  XT(52)=XT(52) + data(t)*sin(MN4 *t);
  XT(53)=XT(53) + data(t)*cos(M4  *t);
  XT(54)=XT(54) + data(t)*sin(M4  *t);
  XT(55)=XT(55) + data(t)*cos(SN4 *t);
  XT(56)=XT(56) + data(t)*sin(SN4 *t);
  XT(57)=XT(57) + data(t)*cos(MS4 *t);
  XT(58)=XT(58) + data(t)*sin(MS4 *t);
end

CoeAB = MT\XT';

% Amplitude
A   = CoeAB(1:2:end-1);
B   = CoeAB(2:2:end);
Amp = sqrt(A.^2 + B.^2);

%  Phase on the beginning of the series
F = atan2(-B,A);

% Phase in degrees
Pha = F.*(180/pi);
Pha(Pha<0) = Pha(Pha <0)+360;

% reconstruct the series
W = w(1:2:end);
for k = 1:nhrs
    Tid(:,k) = Amp(:).*cos(W(:).*k + F(:));
end

%  Tid = M2+S2+N2+K2+K1+O1+P1+Q1;

clc

fprintf('\r\n ');
fprintf('\r\n ');
fprintf('Mean level = '); fprintf('%8.2f', mean(data));
fprintf('\r\n ');
fprintf('Tidal component Amplitude Phase \r\n');
fprintf('M2  '); fprintf('%14.2f', Amp(1) ); fprintf('%14.2f\r\n', Pha(1) );
fprintf('S2  '); fprintf('%14.2f', Amp(2) ); fprintf('%14.2f\r\n', Pha(2) );
fprintf('N2  '); fprintf('%14.2f', Amp(3) ); fprintf('%14.2f\r\n', Pha(3) );
fprintf('K2  '); fprintf('%14.2f', Amp(4) ); fprintf('%14.2f\r\n', Pha(4) );
fprintf('K1  '); fprintf('%14.2f', Amp(5) ); fprintf('%14.2f\r\n', Pha(5) );
fprintf('O1  '); fprintf('%14.2f', Amp(6) ); fprintf('%14.2f\r\n', Pha(6) );
fprintf('P1  '); fprintf('%14.2f', Amp(7) ); fprintf('%14.2f\r\n', Pha(7) );
fprintf('Q1  '); fprintf('%14.2f', Amp(8) ); fprintf('%14.2f\r\n', Pha(8) );