clc;
clear 
%No necessary TOOLBOXES
%% Author: Andres Cornejo, UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO
% Copyright ©2020 Version (09/2020)
%% ****************************** Inputs **********************************
Frequency = 50; %GHz V-band Feeder Uplink
%**************************************************************************

%% ******************* Compute Rain Attenuation, CCDF *********************
ccdf = [kron(10.^(-2:0),[1 2 3 5]),10];
Ar = zeros(1,size(ccdf,2));
for j = 1:1:size(ccdf,2)
    [a,b,~,~,~] = ITU_RP618_12(-92.0,17.113126,-96.74227,1555.0,...
        Frequency,ccdf(1,j));          
    Ar(1,j) = a(1,2);
end
P0 = b;   
% *************************************************************************
sec = 60*60*24*365;      % Total number of seconds per year.
    
%% ********** Rain Attenuation Time-Series by ONERA® Function *************
%Sampling every 60 seconds or 1 minute        
[ARain,Time,lnm,lns]=itur_p1853_0_annex1_s2(ccdf,Ar,P0,60,sec);   
%**************************************************************************