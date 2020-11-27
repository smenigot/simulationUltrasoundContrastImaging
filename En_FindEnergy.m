function [E,bubble,medium,pulse1] = En_FindEnergy(pulse, bubble, medium);
% Computation of the CTR from the propagation result
%
% By Sébastien Ménigot (Université François-Rabelais de Tours)
% sebastien.menigot@univ-tours.fr
%
% sous licence CeCILL (http://www.cecill.info/)


if ~medium.flag
%% Bubble properties
    bubble.name = 'Bubble';
    bubble.centre_bulle = 15e-3;           % [m]    Distance between the probe and the bubbles

    % Propriete de la coque de la bulle
    bubble.ds = 1.0e-9;                     % [m]   Shell thickness (Epaisseur de la coque)
    bubble.Gs = 46e6;                       % [Pa]  Shell shear modulus
    bubble.es = 1;                          % [Pas] Shell shear viscosity
    
    bubble.rho  = 1050;      % [kg/m3]  Density of kidney
    bubble.c    = 1578;      % [m/s]    Speed of sound in kidney
    bubble.p0   = 1.013e5;     % [Pa]     Ambient pressure
    bubble.eL   = 4.0e-3;      % [Pas]    Viscosity in liquid
    bubble.gamma= 660/603;   % [1]      Cp/Cv (SF6)
    
    bubble.ksi          = 0.38;             %       Shell elasticity
    bubble.sigma_break	= 0.13;
    bubble.sigma_liquid = 0.058;            %       
    bubble.ba           = 6.75;             %       B/A nonlinear parameter
    
    bubble.number    = 1;                                    %       number of bubble
    bubble.a0        = 2.5/2 * ones(1,bubble.number) * 1e-6; % [m]   radius
    bubble.R_rupture = 2*bubble.a0;                          % [m]   rupture radius

%% parameter of the first computation
    attn.flag_calcul = true;
    map.flag = true;
else
    attn = medium.attn;
    attn.flag_calcul = false;
    
    map = medium.map;
    map.flag = false;
end

%% Propagation
[time1, echo1, map1, attn1, pulse1] = En_propagation(pulse, bubble, map, attn);

medium.flag = true;
medium.map = map1;
medium.attn = attn1;

pulse.A = -pulse.A;

[time2, echo2, map2, attn2] = En_propagation(pulse, bubble, map, attn1);

%% Power and CTR
dt = time1(2)-time1(1);
% Power without bubble
sum_linear = echo1(:,1) + echo2(:,1);
E.E_linear = sum(sum_linear.^2)/( length(sum_linear)*dt );
% Power with bubble
sum_nonlinear = echo1(:,2) + echo2(:,2);
E.E_nonlinear = sum(sum_nonlinear.^2)/( length(sum_nonlinear)*dt );
% CTR
E.CTR = E.E_nonlinear ./ E.E_linear;

return

