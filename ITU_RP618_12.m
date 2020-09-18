function [output1,output2,output3,output4,output5] = ITU_RP618_12(pos_sat,Latitude,Longitude,Height,Frequency,p)
%No necessary TOOLBOXES
%% Author: Andres Cornejo, UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO
% Copyright ©2020
%% ******* FUNCTION TO EVALUATE RAINFALL AND PERCENTAGE PROBABILITY *******
% Percentage probability of rain in an average year (%)
h0 = importdata('h0.txt');
Lat = importdata('Lat.txt');
Lon = importdata('Lon.txt');

Lat = (Lat(:,1));                   %Latitude ranges -90 a 90 deg.
Lon = -180.+(Lon(1,:));             %Longitude ranges -180 a 180 deg.

for i = 1:1:size(Lat,1)
    if Lat(i,1) < Latitude
    else
        x=i;
    end
end

for j = 1:1:size(Lon,2)
    if Lon(1,j) > Longitude
    else
        y=j;
    end
end

%% ******************** BILINEAR INTERPOLATION ****************************
% ITU-R P.1144-9: Guide to the Application of the Propagation Methods of 
% Radiocommunication Study Group 3 
%**************************************************************************

%% ****************************** STEP 1 **********************************
r = (Latitude);
R = (Lat(x,1));
R_1 = (Lat(x+1,1));
c = (Longitude);
C= (Lon(1,y));
C_1 = (Lon(1,y+1));

if R_1 > R && C_1 > C
    disp('caso 1')
elseif R_1 < R && C_1 > C
%Coordinates Normalization
    n = (r-R_1)/(R-R_1);
    m = (c-C)/(C_1-C);
    R = 0;
    r = n;
    R_1 = 1;
    C = 0;
    c = m;
    C_1 = 1;
elseif R_1 > R && C_1 < C
    disp('caso 3')
elseif R_1 < R && C_1 < C
    disp('caso 4')
end

H0 = h0(x,y)*(R_1-r)*(C_1-c)+h0(x+1,y)*(r-R)*(C_1-c)+h0(x,y+1)*(R_1-r)...
    *(c-C)+h0(x+1,y+1)*(r-R)*(c-C);
hr = H0+0.36; %
%**************************************************************************

%% ****************************** STEP 2 **********************************
estructura = elazr(pos_sat,Latitude,Longitude,35786,6378.14);
if estructura.Elevacion >= 5
    hs = Height/1000;  %Height above the sea level, km.
    Ls = (hr-hs)/sin(estructura.Elevacion*pi/180);
else
    hs = Height/1000;
    Ls = 2*(hr-hs)/((sin(estructura.Elevacion*pi/180)^2)+...
        ((2*(hr-hs)/6378.14)^(1/2))+(sin(estructura.Elevacion*pi/180)));
end
% If hr-hs < 0 para que A = 0 dB ******************************************
if (hr-hs) <= 0
    goto end 
end
%**************************************************************************

%% ****************************** STEP 3 **********************************
% Horizontal projection ***************************************************
LG = Ls*cos(estructura.Elevacion*pi/180);
%**************************************************************************

%% ****************************** STEP 4 **********************************
%Recomendacion ITU-R P.837. Percentage 0.01%
%% ***************** ITU-R P.837-7 ANNEX 1 by ONERA® **********************
[Rp,P0,Rp_m,MT,T] = itur_p0837_7_annex1(0.01,Latitude,Longitude,0);
if Rp == 0
    goto end
end
%**************************************************************************

%% ****************************** STEP 5 **********************************
% Recomendacion ITU-R P.838-3
YR = ITU_RP838(Frequency,estructura.Elevacion,'circular',Rp);
%**************************************************************************

%% ****************************** STEP 6 **********************************
r_001 = 1/(1+(0.78*sqrt(LG*YR/Frequency))-0.38*(1-exp(-2*LG)));
%**************************************************************************

%% ****************************** STEP 7 **********************************
sigma = atan((hr-hs)/(LG*r_001));   %Degress

if sigma > estructura.Elevacion
    LR = LG*r_001/cos(estructura.Elevacion*pi/180);
else
    LR = (hr-hs)/sin(estructura.Elevacion*pi/180);
end

if abs(Latitude) < 36
    X = 36 - abs(Latitude);
else
    X = 0;
end

v_001 = 1/(1+sqrt(sin(estructura.Elevacion*pi/180))*(31*(1-exp...
    (-(estructura.Elevacion/(1+X))))*(sqrt(LR*YR)/(Frequency^2))-0.45));
%**************************************************************************

%% ****************************** STEP 8 **********************************
% Effective Length, LE, of the path:
LE = LR*v_001;
%**************************************************************************

%% ****************************** STEP 9 **********************************
% Attenuation for 0.01% of a year,
% A0.01%,
A_001 = YR*LE;
%**************************************************************************

%% ***************************** STEP 10 **********************************
if p >= 1 || abs(Latitude) >= 36
    B = 0;
elseif p < 1 && abs(Latitude) < 36 && estructura.Elevacion >= 25
    B = -0.005*(abs(Latitude)-36);
else
    B = -0.005*(abs(Latitude)-36)+1.8-4.25*sin(estructura.Elevacion*pi/180);
end

Ap = A_001*(p/0.01)^(-(0.655+0.033*log(p)-0.045*log(A_001)-B*(1-p)...
    *sin(estructura.Elevacion*pi/180)));
%**************************************************************************

%% **************************** RESULTADOS ********************************
output5 = Rp_m; %rainfall rate exceeded for p% of an average month for all months  (mm/h) (12xlength(p))
output4 = MT;   %monthly rain amount from the maps (mm)
output3 = T;    %monthly ground temperature from the maps (Celsius)
output2 = P0;   %percentage probability of rain in an average year (%)
output1 = [A_001,Ap];
%% ************************************************************************
end
