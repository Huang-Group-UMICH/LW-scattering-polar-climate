clc;clear;close all;
addpath(genpath('export_fig'));

%% Loading data

% Path of cloud optical properties
Path_Ice = 'MC6_Ice_OptProp.txt'; % path of Modis collection 6 ice cloud bulk scattering properties
Path_Water = 'MC6_Water_OptProp.txt'; % path of Modis collection 6 water cloud bulk scattering properties

% Load bulk optical properties
OptIce = importdata(Path_Ice, ' ', 1);
OptIce = OptIce.data;
OptWat = importdata(Path_Water, ' ', 1);
OptWat = OptWat.data;

% MC6 ice cloud optical properties
Ice_NDeff = 18; % Number of effective diameter
Ice_NWave = 961; % Number of wavelength
Ice_Deff = zeros(1,Ice_NDeff);
Ice_Wave = zeros(1,Ice_NWave);
Ice_Qext = zeros(Ice_NDeff,Ice_NWave);
Ice_Qabs = zeros(Ice_NDeff,Ice_NWave);
Ice_Qsca = zeros(Ice_NDeff,Ice_NWave);
Ice_Albedo = zeros(Ice_NDeff,Ice_NWave);
Ice_Gfactor = zeros(Ice_NDeff,Ice_NWave);
Ice_Mext = zeros(Ice_NDeff,Ice_NWave);
Ice_Mabs = zeros(Ice_NDeff,Ice_NWave);
for iDeff = 1:Ice_NDeff % Loop over effective diameter from 10 to 180 um with interval 10 um
    for iWave = 1:Ice_NWave % Loop over wavelength (um) from 3.0 to 99.9 um with interval 0.1 um
        TmpIndex = (iDeff-1)*Ice_NWave + iWave;
        Ice_Deff(iDeff)          = OptIce(TmpIndex,1); % Effective diameter (um)
        Ice_Wave(iWave)          = OptIce(TmpIndex,2); % Wavelength (um)
        Ice_Qext(iDeff,iWave)    = OptIce(TmpIndex,3); % Mean extinction efficiency
        Ice_Qabs(iDeff,iWave)    = OptIce(TmpIndex,4); % Mean absorption efficiency
        Ice_Qsca(iDeff,iWave)    = OptIce(TmpIndex,5); % Mean scattering efficiency
        Ice_Albedo(iDeff,iWave)  = OptIce(TmpIndex,6); % Single-scattering albedo
        Ice_Parea                = OptIce(TmpIndex,7); % Projection area (um^2/particle)
        Ice_Volume               = OptIce(TmpIndex,8); % Volume (um^3/particle)
        Ice_Gfactor(iDeff,iWave) = OptIce(TmpIndex,9); % Asymmetry factor
        IWC                      = 0.971 * Ice_Volume; % Ice water content (IWC, 1E-12 g/particle; ice density, 0.917 g/cm^3)
        Ice_Mext(iDeff,iWave)    = Ice_Qext(iDeff,iWave) * Ice_Parea / IWC; % Mass extinction coefficient (m^2/g)
        Ice_Mabs(iDeff,iWave)    = Ice_Qabs(iDeff,iWave) * Ice_Parea / IWC; % Mass absorption coefficient (m^2/g)
    end
end

% MC6 ice cloud optical properties
Wat_NDeff = 18; % Number of effective diameter
Wat_NWave = 509; % Number of wavelength
Wat_Deff = zeros(1,Wat_NDeff);
Wat_Wave = zeros(1,Wat_NWave);
Wat_Qext = zeros(Wat_NDeff,Wat_NWave);
Wat_Qabs = zeros(Wat_NDeff,Wat_NWave);
Wat_Qsca = zeros(Wat_NDeff,Wat_NWave);
Wat_Albedo = zeros(Wat_NDeff,Wat_NWave);
Wat_Gfactor = zeros(Wat_NDeff,Wat_NWave);
Wat_Mext = zeros(Wat_NDeff,Wat_NWave);
Wat_Mabs = zeros(Wat_NDeff,Wat_NWave);
for iDeff = 1:Wat_NDeff % Loop over effective diameter from 10 to 180 um with interval 10 um
    for iWave = 1:Wat_NWave % Loop over wavelength (um) from 3.0 to 99.9 um with interval 0.1 um
        TmpIndex = (iDeff-1)*Wat_NWave + iWave;
        Wat_Deff(iDeff)          = OptWat(TmpIndex,1); % Effective diameter (um)
        Wat_Wave(iWave)          = OptWat(TmpIndex,2); % Wavelength (um)
        Wat_Qext(iDeff,iWave)    = OptWat(TmpIndex,3); % Mean extinction efficiency
        Wat_Qabs(iDeff,iWave)    = OptWat(TmpIndex,4); % Mean absorption efficiency
        Wat_Qsca(iDeff,iWave)    = OptWat(TmpIndex,5); % Mean scattering efficiency
        Wat_Albedo(iDeff,iWave)  = OptWat(TmpIndex,6); % Single-scattering albedo
        Wat_Parea                = OptWat(TmpIndex,7); % Projection area (um^2/particle)
        Wat_Volume               = OptWat(TmpIndex,8); % Volume (um^3/particle)
        Wat_Gfactor(iDeff,iWave) = OptWat(TmpIndex,9); % Asymmetry factor
        LWC                      = 1.0 * Wat_Volume;   % Liquid water content (LWC, 1E-12 g/particle; liquid water density (4 degree), 1.0 g/cm^3)
        Wat_Mext(iDeff,iWave)    = Wat_Qext(iDeff,iWave) * Wat_Parea / LWC; % Mass extinction coefficient (m^2/g)
        Wat_Mabs(iDeff,iWave)    = Wat_Qabs(iDeff,iWave) * Wat_Parea / LWC; % Mass absorption coefficient (m^2/g)
    end
end

%% Plot cloud optical properties

% Figure parameters
ImgShow     = 'on';
ImgSize     = [1536, 1024].*0.8; % Figure size
FSize       = 26;                % Font Size
LegendFSize = 26;                % Legend Font Size
LineWid     = 3;                 % Line width
BoxLineWid  = 3;                 % Box Line width
MakerSize   = 12;                % Marker Size
Marker      = '*';               % Marker

% Figure 1. plot extinction efficiency
TheFig(1) = FigCreat( 1, ImgSize, ImgShow ); % Create new figure for number of Index in size ImgSize
set(TheFig(1), 'Renderer', 'opengl');
plot( 10000./Ice_Wave, Ice_Qext(2,:), 'Color', [0 0 0], 'Linestyle', '-',  'LineWidth', LineWid); hold on  % Ice cloud effective diameter 20 um
plot( 10000./Ice_Wave, Ice_Qext(6,:), 'Color', [0 0 0], 'Linestyle', '-.', 'LineWidth', LineWid);          % Ice cloud effective diameter 60 um
plot( 10000./Wat_Wave, Wat_Qext(2,:), 'Color', [0 0 0], 'Linestyle', ':',  'LineWidth', LineWid); hold off % Water cloud effective diameter 20 um
xlim([0, 3250]);
ylim([0.5, 3]);
xlabel('Wavenumber (cm^{-1})','FontSize',FSize)
ylabel('Extinction efficiency','FontSize',FSize)
set(gca,'FontSize', FSize, ...
        'XTick', 0:500:3500, ...
        'XTicklabel', 0:500:3500, ...
        'YTick', 0.5:0.5:3, ...
        'YTicklabel', 0.5:0.5:3, ...
        'TickLength', 2.5*get(gca,'ticklength'), ...
        'Linewidth', BoxLineWid, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on');
Handle = legend( ['Ice(20',char(181),'m)'], ['Ice(60',char(181),'m)'], ['Water(20',char(181),'m)']);
set( Handle, 'FontSize', LegendFSize, 'Box', 'off', 'Color', 'none', 'Position', [0.50, 0.25, .25, .25]);

% Figure 2. plot single-scattering albedo
TheFig(2) = FigCreat( 2, ImgSize, ImgShow ); % Create new figure for number of Index in size ImgSize
set(TheFig(2), 'Renderer', 'opengl');
plot( 10000./Ice_Wave, Ice_Albedo(2,:), 'Color', [0 0 0], 'Linestyle', '-',  'LineWidth', LineWid); hold on  % Ice cloud effective diameter 20 um
plot( 10000./Ice_Wave, Ice_Albedo(6,:), 'Color', [0 0 0], 'Linestyle', '-.', 'LineWidth', LineWid);          % Ice cloud effective diameter 60 um
plot( 10000./Wat_Wave, Wat_Albedo(2,:), 'Color', [0 0 0], 'Linestyle', ':',  'LineWidth', LineWid); hold off % Water cloud effective diameter 20 um
xlim([0, 3250]);
ylim([0.2, 1]);
xlabel('Wavenumber (cm^{-1})','FontSize',FSize)
ylabel('Single-scattering albedo','FontSize',FSize)
set(gca,'FontSize', FSize, ...
        'XTick', 0:500:3500, ...
        'XTicklabel', 0:500:3500, ...
        'YTick', 0.2:0.2:1, ...
        'YTicklabel', 0.2:0.2:1, ...
        'TickLength', 2.5*get(gca,'ticklength'), ...
        'Linewidth', BoxLineWid, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on');
Handle = legend( ['Ice(20',char(181),'m)'], ['Ice(60',char(181),'m)'], ['Water(20',char(181),'m)']);
set( Handle, 'FontSize', LegendFSize, 'Box', 'off', 'Color', 'none', 'Position', [0.60, 0.18, .25, .25]);

% Figure 3. plot asymmetry factor
TheFig(3) = FigCreat( 3, ImgSize, ImgShow ); % Create new figure for number of Index in size ImgSize
set(TheFig(3), 'Renderer', 'opengl');
plot( 10000./Ice_Wave, Ice_Gfactor(2,:), 'Color', [0 0 0], 'Linestyle', '-',  'LineWidth', LineWid); hold on  % Ice cloud effective diameter 20 um
plot( 10000./Ice_Wave, Ice_Gfactor(6,:), 'Color', [0 0 0], 'Linestyle', '-.', 'LineWidth', LineWid);          % Ice cloud effective diameter 60 um
plot( 10000./Wat_Wave, Wat_Gfactor(2,:), 'Color', [0 0 0], 'Linestyle', ':',  'LineWidth', LineWid); hold off % Water cloud effective diameter 20 um
xlim([0, 3250]);
ylim([0.2, 1]);
xlabel('Wavenumber (cm^{-1})','FontSize',FSize)
ylabel('Asymmetry factor','FontSize',FSize)
set(gca,'FontSize', FSize, ...
        'XTick', 0:500:3500, ...
        'XTicklabel', 0:500:3500, ...
        'YTick', 0.2:0.2:1, ...
        'YTicklabel', 0.2:0.2:1, ...
        'TickLength', 2.5*get(gca,'ticklength'), ...
        'Linewidth', BoxLineWid, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on');
Handle = legend( ['Ice(20',char(181),'m)'], ['Ice(60',char(181),'m)'], ['Water(20',char(181),'m)']);
set( Handle, 'FontSize', LegendFSize, 'Box', 'off', 'Color', 'none', 'Position', [0.50, 0.30, .25, .25]);

% Save figures
export_fig(TheFig(1),'Figures/ExtinctionEfficiency','-pdf','-r150','-q101','-a4','-opengl');
export_fig(TheFig(1),'Figures/ExtinctionEfficiency','-png','-tiff','-r150','-q101','-a4','-opengl');
export_fig(TheFig(2),'Figures/SingleScatteringAlbedo','-pdf','-r150','-q101','-a4','-opengl');
export_fig(TheFig(2),'Figures/SingleScatteringAlbedo','-png','-tiff','-r150','-q101','-a4','-opengl');
export_fig(TheFig(3),'Figures/AsymmetryFactor','-pdf','-r150','-q101','-a4','-opengl');
export_fig(TheFig(3),'Figures/AsymmetryFactor','-png','-tiff','-r150','-q101','-a4','-opengl');

