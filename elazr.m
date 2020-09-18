function outputstruct = elazr(pos_sat,latitude,longitude,RO,Re)
%No necessary TOOLBOXES
%% Author: Andres Cornejo, UNIVERSIDAD NACIONAL AUTONOMA DE MEXICO
% Copyright ©2020
%% ************************************************************************
% Function to determine Elevation, Azimuth and Range
% Satellite Position, GEO = pos_sat [deg]
% GW Latitude = latitude [deg]
% GW Longitude = longitude [deg]
% RO = 35786 [Km]; Range Earth-satellite
% Re = 6378.14 [Km]; Radious of the Earth
% *************************************************************************
%% ************************************************************************
narginchk(5, 5);  
phi_rad = latitude*pi/180;
delta_rad = (longitude*pi/180)-(pos_sat*pi/180);
Cos_PHI = cos(phi_rad)*cos(delta_rad);

R = 42643.7*sqrt(1-0.29577*Cos_PHI);
E = atan((Cos_PHI-(Re/(Re+RO)))/sqrt(1-Cos_PHI^2))*(180/pi);

Az = azimuth(latitude,longitude,0,pos_sat,'deg');

outputstruct = struct('R',R,'Elevacion',E,'Azimuth',Az);

end
%**************************************************************************