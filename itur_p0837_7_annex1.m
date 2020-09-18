function [Rp,P0,Rp_m,MT,T] = itur_p0837_7_annex1(p,lat,lon,bool_plot,model_maps)
%
% function [Rp,P0,Rp_m,MT,T] = itur_p0837_7_annex1(p,lat,lon,bool_plot,model_maps)
%
% This function calculates the rainfall rate exceeded for p% of an average
% year or month
%
% Input parameters :
%   p           = probability level (%)
%   lat         = latitude          [+90 .. -90] (deg)
%   lon         = longitude         [-180E .. 180E] (deg)
%   bool_plot   = Boolean if user wants to compute and plot monthly results (0:No - 1:Yes)
%   model_maps  = if the digital maps are already loaded into the
%                 workspace, the user can save time by using them
%                 directly (if this input parameter is empty, the software will automatically load the data)
%
% Output parameters :
%   Rp      = rainfall rate exceeded for p% of an average year  (mm/h)
%   p0      = percentage probability of rain in an average year (%)
%   Rp_m    = rainfall rate exceeded for p% of an average month for all months  (mm/h) (12xlength(p))
%   MT      = monthly rain amount from the maps (mm)
%   T       = monthly ground temperature from the maps (°C)
%
% Note :
%   No
%
% Called functions :
%   No
%
% Necessary toolboxes :
%   No
%
% by Xavier Boulanger
% ONERA, France
% Any questions : email <Xavier.Boulanger@onera.fr>
%
% rel. 0.0
% release history:
%   0.0 (09/2016)    - original version
%-------------------------------------------------------------------------------------------------------------------

if lon>180
    lon = lon - 360;
end

if nargin<5
    disp('Loading maps in progress...')
    model_maps = load('_ITU-R_P0837-7_MT_T_maps.mat');
    disp('Maps loaded')
end
R001 = load('_ITU-R_P0837-7_R001_map.mat');

sig         = 1.26;
N_ii        = [31 28.25 31 30 31 30 31 31 30 31 30 31];

for ii=1:12
    MT_ii(ii)  = interp2(model_maps.lon_025,model_maps.lat_025,model_maps.Mt_025(:,:,ii),lon,lat,'bilinear');   % Montly precipitation amount
    T_ii(ii)   = interp2(model_maps.lon_075,model_maps.lat_075,model_maps.T_075(:,:,ii),lon,lat,'bilinear')-273.15;    % Monthly temperature

    r_ii(ii)   = 0.5874*exp(0.0883*T_ii(ii));
    if T_ii(ii) < 0; r_ii(ii) = 0.5874; end

    P0_ii(ii)  = 100*MT_ii(ii)/(24*N_ii(ii)*r_ii(ii));                               % Monthly probability of rain
    if P0_ii(ii)>70;
        P0_ii(ii) = 70;
        r_ii(ii) = (100*MT_ii(ii))./(70*24*N_ii(ii));
    end
end
P0 = sum(P0_ii.*N_ii)/sum(N_ii); % Yearly CCDF of rainfall rate

for ip = 1:length(p)
    if p(ip)>P0
        Rp(ip)=0;
    else
        R_low   = 0;
        R_high  = 500;
        R_ref   = interp2(R001.longitude,R001.latitude,R001.R001,lon,lat);if R_ref == 0, R_ref = 10; end

        for ii=1:12
            Pc_ii(ii) = 0.5*erfc((log(R_ref)+sig^2/2-log(r_ii(ii)))/(sig*sqrt(2)));  % Monthly conditional CCDF of rainfall rate
            P_ii(ii)  = P0_ii(ii)*Pc_ii(ii);
        end
        P_fin   = sum(P_ii.*N_ii)/sum(N_ii); % Yearly CCDF of rainfall rate

        rel_err = 100*(P_fin-p(ip))./p(ip);

        while abs(rel_err)>1e-3 % (in %)
            if P_fin<p(ip)
                R_high = R_ref;
            else
                R_low  = R_ref;
            end
            R_ref  = (R_high + R_low)./2;

            for ii=1:12
                Pc_ii(ii) = 0.5*erfc((log(R_ref)+sig^2/2-log(r_ii(ii)))/(sig*sqrt(2)));  % Monthly conditional CCDF of rainfall rate
                P_ii(ii)  = P0_ii(ii)*Pc_ii(ii);
            end
            P_fin   = sum(P_ii.*N_ii)/sum(N_ii); % Yearly CCDF of rainfall rate

            rel_err = 100*(P_fin-p(ip))./p(ip);
        end
        Rp(ip) = R_ref;
    end
end

MT = MT_ii;
T  = T_ii;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTATION OF MONTHLY RESULTS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bool_plot
    ind=0;
    for ii=1:12
        for ip = 1:length(p)
            ind=ind+1;
            clc;disp(['Computation of monthly results in progress: ' num2str(ind./(12*length(p))*100) ' %']);
            if p(ip)>P0_ii(ii)
                Rp_m(ii,ip)=0;
            else
                R_low   = 0;
                R_high  = 500;
                R_ref   = interp2(R001.longitude,R001.latitude,R001.R001,lon,lat);if R_ref == 0, R_ref = 10; end

                Pc_ii = 0.5*erfc((log(R_ref)+sig^2/2-log(r_ii(ii)))/(sig*sqrt(2)));  % Monthly conditional CCDF of rainfall rate
                P_ii  = P0_ii(ii)*Pc_ii;

                rel_err = 100*(P_ii-p(ip))./p(ip);

                while abs(rel_err)>1e-3 % (in %)
                    if P_ii<p(ip)
                        R_high = R_ref;
                    else
                        R_low  = R_ref;
                    end
                    R_ref  = (R_high + R_low)./2;

                    Pc_ii = 0.5*erfc((log(R_ref)+sig^2/2-log(r_ii(ii)))/(sig*sqrt(2)));  % Monthly conditional CCDF of rainfall rate
                    P_ii  = P0_ii(ii)*Pc_ii;

                    rel_err = 100*(P_ii-p(ip))./p(ip);
                end
                Rp_m(ii,ip) = R_ref;
            end
        end
    end

    months={'January','February','March','April'...
        'May','June','July','August',...
        'September','October','November','December','Total'};
    %% Figures
    h=figure('Position',[0 0 1600 800],'PaperPositionMode','auto');
    subplot(1,2,1)
    set(gca,'fontsize',24);
    bar([1:12],MT_ii); grid on
    xlabel('Months'); ylabel('Rain amount (mm)'); xlim([0 13])
    subplot(1,2,2)
    set(gca,'fontsize',24)
    plot([1:12],T_ii,'r-o'); grid on
    xlabel('Months'); ylabel('Temperature (°C)'); xlim([0 13])

    color_plot = {'k','k--','k-.','g','g--','g-.','r','r--','r-.','b','b--','b-.'};
    h=figure('Position',[0 0 1200 800],'PaperPositionMode','auto');
    set(gca,'fontsize',20);
    for ii=1:12
        semilogx(p,Rp_m(ii,:),color_plot{ii},'linewidth',2); hold on
    end
    semilogx(p,Rp,'k','linewidth',3);hold on
    xlabel('P(R>R*) (%)'); ylabel('Rainfall rate exceeded R* (mm/h)'); grid on
    xlim([1e-2 50])
    legend(months)
else
    Rp_m = 'Not computed';
end

return