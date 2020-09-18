function [Arain, Time,lns,lnm]=itur_p1853_0_annex1_s2(ccdf,Ar,P0,sampling_rate,duration,seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function itur_p1853_0_annex1_s2(ccdf,Ar,P0,sampling_rate,duration,seed)
%
% The function itur_p1853_0_annex1_s2 generates rain attenuation time series 
% according to ITU-R recommendation P.1853-0 annex 1 section 2, starting
% from Rain attenuation statistics.
%
% Rain attenuation statistics can be determined from local measured data, or,
% in the absence of measured data, the rain attenuation prediction methods
% in Rec. ITU-R P.530 for terrestrial paths
% or Rec. ITU-R P.618 for Earth-space paths can be used.
% 
% Input parameters :
% ****************
% ccdf :            input CCDF, a suggested set of time percentages is  
%          [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10] (%)
% Ar :              Attenuation levels exceeded for percentages of time 
%                   stored in ccdf                                      (dB)
% P0 :              Probability of the rain on the path. P0 can be well 
%                   approximated as P0(Lat,Lon) deriverd in ITU-R P.837 (%)
% sampling_rate :   Sampling rate                                       (Hz)
% duration :        Duration of the generated time series               (s)
% seed :            Seed to initialize the random generator             (-)
%                   (default value : no specific initialization)
%
% Output parameters :
% *****************
% Arain :           Time series of rain attenuation                     (dB)
% Time :            Time vector corresponding to the time series        (s)
%
% Called functions :
% ****************
% No
%
% Necessary toolboxes :
% *******************
% No
%
% Example :
% *******
% ccdf = [kron(10.^(-2:0),[1 2 3 5]),10];
% Ar = [16.0243, 12.1351, 10.1633, 8.0040, 5.6312, 3.8381, 3.0224, 2.2024, ...
%       1.3946, 0.8555, 0.6334, 0.4271, 0.2434];
% P0 = 4.9095;
% sampling_rate = 1;
% duration = 86400;
% [Arain, Time]=itur_p1853_0_annex1_s2(ccdf,Ar,P0,sampling_rate,duration);
%
%
% by Laurent CASTANET, Guillaume CARRIE & Nicolas Jeannin, ONERA, France
% Any questions : email <Laurent.Castanet@onera.fr>,
% <Guillaume.Carrie@onera.fr>, <Nicolas.jeannin@onera.fr>
% ONERA, October 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Estimation of  the Log-normal and Aoffset parameters of the ONERA-CNES channel model :
% ------------------------------------------------------------------------------------
[lnm,lns,Aoffset] = itur_p1853_0_annex1_s2_2_A_C(ccdf,Ar,P0);


% Execution of the ONERA-CNES rain attenuation synthesiser :
% --------------------------------------------------------
if exist('seed','var') & ~isempty(seed)
   randn('state',seed);
end
beta = 2e-4; % (s-1) value recommended by ITU-R P.1853-0 annex 1 section 2.2 B.
[Arain, Time] = itur_p1853_0_annex1_s2_2_D(sampling_rate,duration,lnm,lns,beta,Aoffset) ;


% Graphical output :
% ----------------
% CCDF Fit
figure
ccdfth = 100*0.5*(1-erf((log(Ar)-lnm)/(lns*sqrt(2))));
semilogx(ccdf,Ar,ccdfth,Ar,'LineWidth',2)
title(['Original and approximated CCDF']);
ylabel('Attenuation Exceeded (dB)')
xlabel('Time (%)')
legend('Original','LogN')

% Synthesized Time-Series
figure ;
plot(Time,Arain,'r-') ; grid on ;
title('\fontsize{14}\bfITU-R P.1853 - Anexo 1') ;
xlabel('\fontsize{14}Time (s)');
ylabel('\fontsize{14}Rain attenuation (dB)');


% -------------------------------------------------------------------------
% *************************     Sub Functions     *************************
% -------------------------------------------------------------------------

function [lnm,lns,Aoffset] = itur_p1853_0_annex1_s2_2_A_C(ccdf,Aexc,P0)
%
% function [lnm,lns,Aoffset] = itur_p1853_0_annex1_s2_2_A_C(ccdf,Aexc,P0)
%
% The function itur_p1853_0_annex1_s2_2_A_C assess the parameters
% m, sigma and Aoffset from a given CCDF of rain attenuation according to
% ITU-R recommendation P.1853-0 annex 1 section 2.2 A & C.
%
% Input
% -----
% ccdf = input CCDF (%)
% Aexc = attenuation levels exceeded for percentages of time stored in ccdf
% (dB)
% P0 = probability to have rain attenuation on the link (%)
%
% Output
% ------
% lnm = m parameter 
% lns = sigma parameter 
% Aoffset = Aoffset parameter 
%

ccdf = ccdf(:);
Aexc = Aexc(:);
ccdf = ccdf(ccdf<=100);
Aexc = Aexc(ccdf<=100);

% Log-normal fit of the input CCDF between pinf and psup
ind_perc = (ccdf>=0 & ccdf<=P0);
coeff = polyfit(log(Aexc(ind_perc)),-sqrt(2).*erfinv(1-2*ccdf(ind_perc)/100),1);
lnm = -coeff(2)/coeff(1);
lns = -1/coeff(1);

% Assessment of the attenuation offset Aoffset correponding to P0 in the input CCDF
Aoffset = exp(lnm + lns*sqrt(2)*erfinv(1-2*P0/100));
% -------------------------------------------------------------------------


function [att,time] = itur_p1853_0_annex1_s2_2_D(Ts,D,lnm,lns,beta,Aoffset)
%
% function [att,time] = itur_p1853_0_annex1_s2_2_D(Ts,D,beta,lnm,lns,Aoffset)
%
% The function itur_p1853_0_annex1_s2_2_D synthesises rain attenuation time series
% according to ITU-R recommendation P.1853-0 annex 1 section 2.2 D.
%
% Input
% -----
%
% Ts :              Sampling time                                           (s)
% D :               Duration of the time series to be synthesised           (s)
% lnm :             m parameter                                             (-)
% lns :             sigma parameter                                         (-)
% beta :            Beta parameter (default value : 2e-4)                   (s-1)
% Aoffset :         Aoffset parameter                                       (dB)
%
% Output
% ------
%
% att :             Rain attenuation time series                            (dB, positive values)
% time :            Corresponding time vector                               (s)
%

% Time added at the beginning of the synthesised WGN (then removed from the
% final time series)
Nadd = 200000;  

% Number of samples to synthesise
Ns = floor(D/Ts)+1;

n1 = randn(1,Ns+Nadd);		% Random synthesis of the WGN

% Definition of the low-pass filter
ro = exp(-beta*Ts);
b = sqrt(1-ro^2);
a = [1 -ro]; 

% Low-pass filtering of the WGN
Y = filter(b,a,n1); 
Y(1:Nadd) = [];

% Log-normal ponderation of the WGN
att = exp(lnm+lns.*Y);

% Subtraction of the attenuation offset to the time series
att = max(att - Aoffset,0);

time = (0:Ns-1)*Ts;
% -------------------------------------------------------------------------
