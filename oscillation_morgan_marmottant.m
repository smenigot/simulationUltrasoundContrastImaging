% Equation of Marmottant (oscillation)
%
%   [1] Marmottant, P.; van der Meer, S.; Emmer, M.; Versluis, M.; de Jong, N.; Hilgenfeldt, S. & Lohse, D.
%       A model for large amplitude oscillations of coated bubbles accounting for buckling and rupture
%       Journal of the Acoustical Society of America, 2005, 118, 3499-3505
%
% By Sébastien Ménigot (Université François-Rabelais de Tours)
% sebastien.menigot@univ-tours.fr
%
% sous licence CeCILL (http://www.cecill.info/)

function dy=oscillation_morgan_marmottant(t,y,pulse,bubble,destruction)

temps = pulse.temps';
pulsation = pulse.pression';
R_buckling = bubble.R_buckling;
R0 = bubble.a0;
ksi = bubble.ksi;
rho = bubble.rho;
P0 = bubble.p0;
eL = bubble.eL;
es = bubble.es;
ds = bubble.ds;
gamma = bubble.gamma;
c = bubble.c;
sigma_liquid = bubble.sigma_liquid;
R_break = bubble.R_break;
%%
dy(1)=y(2);
val=interp1(temps,pulsation,t','linear');
%%
if ~destruction
    if y(1) <=R_buckling/R0
        sigma_R=0;
    elseif y(1)>R_buckling/R0
        sigma_R=ksi*((abs(y(1))*R0/R_buckling).^2-1);
    end
else
    if y(1) <=R_buckling/R0
        sigma_R=0;
    elseif y(1)>R_buckling/R0  && (ksi*((y(1)*R0/R_buckling).^2-1))<sigma_liquid
        sigma_R=ksi*((y(1)*R0/R_buckling).^2-1);
    elseif y(1)>R_break/R0  || (ksi*((y(1)*R0/R_buckling).^2-1))>=sigma_liquid
        sigma_R=sigma_liquid;
    end
end

%% Marmottant
sigma_R0=ksi*((R0/R_buckling)^2-1);
kappa_s=2.4e-9;

dy(2)=abs(y(1))^(-1) * ...
        ( -3/2*y(2)*y(2) + ...
            1/(rho*R0*R0)* ( ...
                ( P0 + 2*sigma_R0/R0 ) * abs(y(1))^(-3*gamma) * ( 1-3*gamma*R0/c*y(2) ) - ...
                    2*sigma_R/abs(y(1)) * (1/R0) - ...
                    4*y(2)/abs(y(1)) * ( eL + kappa_s / (R0*abs(y(1))) ) - ...
                    (P0+val)...
                 )...
         );
%%
dy=dy';

