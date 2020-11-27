% Propabulle est un programme qui calcul l'�nergie r�trodiffus�e par le
% tissu et le nuage de microbulles dilu� � 1/2000.
%
% By S�bastien M�nigot (Universit� Fran�ois-Rabelais de Tours)
% sebastien.menigot@univ-tours.fr
%
% sous licence CeCILL (http://www.cecill.info/)

close all;
clear;
clc;
warning off;
%% Sine wave excitation
pulse.Nc = 2.3;   %       Cycle number
pulse.A  = 400e3; % [Pa]  Pulse amplitude
pulse.f0 = 2e6;   % [Hz]  Pulse center frequency

%%
medium.flag = false;    % Medium property : computation ? false : no  -  true : yes
bubble.flag = false;	% Bubble property : computation ? false : no  -  true : yes

E = En_FindEnergy(pulse, bubble, medium);

fprintf('CTR = %2.2f dB\n',10*log10(E.CTR));
