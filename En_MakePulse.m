% function pulse = En_MakePulse( pulse )
% Fonction qui créé l'impulsion
%
% By Sébastien Ménigot (Université François-Rabelais de Tours)
% sebastien.menigot@univ-tours.fr
%
% sous licence CeCILL (http://www.cecill.info/)

function pulse = En_MakePulse( pulse )
warning off;

%% Parametres
A0          = pulse.A;
f0          = pulse.f0;
fs          = pulse.fs;
Nc          = pulse.Nc;

%% --- Time ---
dt = 1/fs;                % Sample interval
tp = pulse.ts;

%% Pulse Ref
fRef = 3e6;
xRef = A0 .* exp(-( (pi*fRef*tp/Nc).^2) ) .*sin(2*pi*fRef.*tp);
PRef = sum(xRef.*xRef)/(length(xRef)/fs);

Gaussienne = exp(-( (pi*fRef*tp/Nc).^2) );

%%
x = A0 .* Gaussienne .*sin(2*pi*f0.*tp);

%% Transducteur
signal_filtre = En_Transducteur(x,fs,'out2');
%% -- Pulse final ---
pulse.t = (0:length(signal_filtre)-1) * dt;
pulse.p = signal_filtre;

return