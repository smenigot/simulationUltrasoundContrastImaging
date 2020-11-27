% function [temps, surface_att] = En_Oscillate(bubble, time_tissue1, echo_tissue1)
% Computation of microbubble oscillation
%
% By Sébastien Ménigot (Université François-Rabelais de Tours)
% sebastien.menigot@univ-tours.fr
%
% sous licence CeCILL (http://www.cecill.info/)

function [temps, surface_att] = En_Oscillate(bubble, time_tissue1, echo_tissue1)

%% Pulse
pulsation.temps = time_tissue1;
pulsation.pression = echo_tissue1;

cond_init1=[1 0];   %   normalized radius and speed
destruction = false;

bubble.R_buckling	= bubble.a0;	%       radius for buckling
bubble.R_break      = bubble.R_buckling * sqrt(1 + bubble.sigma_break / bubble.ksi);
bubble.sigma_R0     = bubble.ksi * ( (bubble.a0/bubble.R_buckling)^2 - 1);

[t1,y1]=ode45(@oscillation_morgan_marmottant,pulsation.temps,cond_init1,[],pulsation,bubble,destruction);

%% Broken case
detruit = find( (y1(:,1).*bubble.a0) > bubble.R_break );
if isempty(detruit) == 0
    destruction = true;
    
    % pulse after broken
    pulse2 = pulsation;
    pulse2.pression = pulsation.pression(detruit(1):end);
    pulse2.temps = pulsation.temps(detruit(1):end);
    
    % initial conditions
    cond_init2=[y1(detruit(1),1) y1(detruit(1),2)];
    
    % resolution of differential equation
    [t2,y2]=ode45(@oscillation_morgan_marmottant,pulse2.temps,cond_init2,[],pulse2,bubble,destruction);
    
    y = bubble.a0 .* [ y1(1:detruit(1)-1,1);y2(1:end,1)]; % rayon de la bulle
    v = bubble.a0 .* [ y1(1:detruit(1)-1,2);y2(1:end,2)]; % vitesse du rayon de la bulle
    
else
    y = bubble.a0 .* y1(:,1);
    v = bubble.a0 .* y1(:,2);
end
%% acceleration
acc = [diff(v)./(t1(2)-t1(1));0];

% Pressure surace
temps = t1;
surface = bubble.rho .* (bubble.centre_bulle).^(-1) .* (y.^2.*acc+2*y.*v.^2);
%% Attenuation of contrast agent
% [1] Gorce, J. M.; Arditi, M. & Schneider, M. 
%     Influence of Bubble Size Distribution on the Echogenicity of 
%     Ultrasound Contrast Agents: A Study of SonoVue. 
%     Investigative radiology, 2000, 35, 661-671
f   =  [0:0.01:0.25              1 1.5 1.625 2   3   4   5   6   7   8   9   10  15  17.5 20:80];
att =  -[zeros(size(0:0.01:0.25)) 3.5 4.4 4.2   4 3.4 2.6 1.9 1.4 1.2 1   0.9 0.8 0.4 0.2  zeros(size(20:80))];
pol = polyfit(f, att, 20);
f2     = (0:ceil(length(surface)/2)-1)/ceil(length(surface)/2)/(temps(2)-temps(1))/1e6;
att_po = polyval(pol,f2);
attenuation_bulle = [att_po att_po(end:-1:1)]';

surface_att = real(ifft( fft(surface) .* 10.^(attenuation_bulle/10))) .* exp(-10^(0.45/10)*temps*bubble.c);









