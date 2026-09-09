% Analysis of Rosario Strait ADCP data collected Oct 2024-Jan 2025
% using Sea Spider seafloor tripod with uplooking Sig250
% Velocity and power statistics for OPALCO tidal project
% J. Thomson, Feb 2025

clear all, close all

scriptdir = fileparts(mfilename('fullpath'));
datadir = fullfile(scriptdir, 'data')

load([datadir '/avgd/S103607A002_SeaSpider_avgd.mat'])

%% turbine params
cutin = 0.7; % m/s
rated = 2.25; % m/s
efficiency = 0.39;
R = 13.5; % rotor radius (m)
hubdepth = 3.5 + R; % hubheight (m) below surface

%% cell height above seabed

SeaSpiderHeight = 0.7; % m
z = SeaSpiderHeight + Config.Average_BlankingDistance + Config.Average_CellSize/2 + Config.Average_CellSize * [1:Config.Average_NCells]; % m above seabed for each cell
z = single(z); % change from unit64 to single precision, for later opperations

%% rename variables (for convenience)

t = Data.Average_Time; % matlab datenum (days from 1-1-0000)
u = Data.Average_VelEast; % m/s eastward
v = Data.Average_VelNorth; % m/s northward
w = ( Data.Average_VelUp1 + Data.Average_VelUp2) ./ 2; % up
waterdepth = Data.Average_Pressure + SeaSpiderHeight; % m total water depth

%% remove data when out of water

ondeck = find( waterdepth < 10 );
u(ondeck,:) = NaN;
v(ondeck,:) = NaN;
w(ondeck,:) = NaN;

%% remove surface reflections

maxrange = (waterdepth-SeaSpiderHeight) * cosd(20); % geometric max range, very conservative, lose the upper ~6 m
for i=1:length(t)
    maxbin(i) = find( (z-SeaSpiderHeight) > maxrange(i),1);
    u(i, maxbin(i):end) = NaN;
    v(i, maxbin(i):end) = NaN;
    w(i, maxbin(i):end) = NaN;
end


%% look at quality metrics / ancillary data
inwater = find( waterdepth > 20 );

% orientation
fh1 = figure(1); clf
subplot(3,1,1), plot(t(inwater),Data.Average_Pitch(inwater), t(inwater),Data.Average_Roll(inwater)), datetick, legend('pitch','roll'), set(gca,'YLim',[-10 10])
subplot(3,1,2), plot(t(inwater),Data.Average_Heading(inwater)), datetick, legend('heading'), set(gca,'YLim',[0 360])
subplot(3,1,3), plot(t(inwater),Data.Average_Temperature(inwater)), datetick, legend('Temperature')
safeExport(fh1, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_PRH.png'))

% amplitude
fh2 = figure(2); clf
for i=1:4
    subplot(4,1,i), pcolor(t(inwater), z, eval(['Data.Average_AmpBeam' num2str(i) '(inwater,:)'])'), shading flat, datetick, hold on, axis tight
    plot(t, z(maxbin),'r.', 'linewidth', 1)
    plot(t, waterdepth,'k-', 'linewidth', 1), ylabel('z [m]')
    cb = colorbar; cb.Label.String = 'dB';
end
safeExport(fh2, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_Amp.png'))

% correlation
fh3 = figure(3); clf
for i=1:4
    subplot(4,1,i), pcolor(t(inwater), z, eval(['Data.Average_CorBeam' num2str(i) '(inwater,:)'])'), shading flat, datetick, hold on, axis tight
    plot(t, z(maxbin),'r.', 'linewidth', 1)
    plot(t, waterdepth,'k-', 'linewidth', 1), ylabel('z [m]')
    cb = colorbar; cb.Label.String = '%';
end
safeExport(fh3, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_Cor.png'))

%  difference in vertical velocity estimates
fh4 = figure(4); clf
pcolor(t(inwater), z, (Data.Average_VelUp1(inwater,:) - Data.Average_VelUp2(inwater, :))'), shading flat, datetick, hold on, axis tight
plot(t, z(maxbin),'r.', 'linewidth', 1)
plot(t, waterdepth,'k-', 'linewidth', 1), ylabel('z [m]')
caxis([-0.5 0.5]), cb = colorbar; cb.Label.String = '\Delta w [m/s]';
safeExport(fh4, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_Deltaw.png'))


%% component and speed plots (pcolors)

speed = (u.^2 + v.^2).^.5;

fh5 = figure(5); clf
%colormap(nawhimar)
subplot(3,1,1), pcolor(t, z, u'), shading flat, datetick, hold on, axis tight
caxis([-2 2]), cb = colorbar; cb.Label.String = 'east [m/s]';
plot(t, waterdepth,'k-', 'linewidth', 1), %ylabel('height above seabed, z [m]')
subplot(3,1,2), pcolor(t, z, v'), shading flat, datetick, hold on, axis tight
caxis([-2 2]), cb = colorbar; cb.Label.String = 'north [m/s]';
plot(t, waterdepth,'k-', 'linewidth', 1), ylabel('height above seabed, z [m]')
subplot(3,1,3), pcolor(t, z, w'), shading flat, datetick, hold on, axis tight
caxis([-2 2]), cb = colorbar; cb.Label.String = 'up [m/s]';
plot(t, waterdepth,'k-', 'linewidth', 1), %ylabel('height above seabed, z [m]')
safeExport(fh5, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_UVW.png'))

fh6 = figure(6); clf
colormap(jet)
pcolor(t, z, speed'), shading flat, datetick, hold on, axis tight
caxis([0 3.5]), cb = colorbar; cb.Label.String = 'horizontal speed [m/s]';
plot(t, waterdepth,'k-', 'linewidth', 1), ylabel('z [m]')
safeExport(fh6, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_speed.png'))


%% tidal ellipses and u,v scatter
principalaxis = nanmean( atand( v ./ u ) );
allz = ones(length(t),1)*z;

fh7 = figure(7); clf
colormap copper
scatter(u(:),v(:), 1, allz(:),'filled')
cb = colorbar; cb.Label.String = 'z [m]';
xlabel('u, east [m/s]'), ylabel('v, north [m/s]')
axis([-3 3 -3 3]), grid
safeExport(fh7, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_UVscatter.png'))


%% vertical profiles

fh8 = figure(8); clf
subplot(1,3,1), plot( nanmean(speed), z), hold on, area([0 1.1],[94 94],83), plot([0 1.1],mean(waterdepth(inwater) - hubdepth)*[1 1],'k--'), 
plot([0 1.1],mean(waterdepth(inwater) - hubdepth+R)*[1 1],'k:'), plot([0 1.1],mean(waterdepth(inwater) - hubdepth-R)*[1 1],'k:'), 
ylabel('z [m]'), xlabel('avg speed [m/s]'), 
subplot(1,3,2), plot( 0.5 * 1030 * nanmean( speed.^3 ), z), hold on, area([0 1500],[94 94],83), plot([0 1500],mean(waterdepth(inwater) - hubdepth)*[1 1],'k--'), 
plot([0 1500],mean(waterdepth(inwater) - hubdepth+R)*[1 1],'k:'), plot([0 1500],mean(waterdepth(inwater) - hubdepth-R)*[1 1],'k:'), 
ylabel('z [m]'), xlabel('avg avail. power density [W/m^2]'), 
subplot(1,3,3), plot( 360+principalaxis, z), hold on, area([300 320],[94 94],83), plot([300 320],mean(waterdepth(inwater) - hubdepth)*[1 1],'k--'), 
plot([300 320],mean(waterdepth(inwater) - hubdepth+R)*[1 1],'k:'), plot([300 320],mean(waterdepth(inwater) - hubdepth-R)*[1 1],'k:'), 
ylabel('z [m]'), xlabel('principal axis [deg M]')
safeExport(fh8, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_profiles.png'))


%% AEP estimate for single O2 turbine, based on turbine parameters at top of script, 
% note that actual production would be double

% hubheight speed
for i=1:length(t)
    [dz zi] = min( abs(z - (waterdepth(i) - hubdepth) ) ) ;
    hubbin(i) = zi;
    hubspeed(i) = speed(i,zi);
end

% annual definition
Tyear = 8760; % hours in a year
dt = median(diff(t)) * 24; % hourly spacing of data
T = range(t(inwater)) * 24; % hours of valid data (103 days)

% power regions
region1 = find( hubspeed < cutin );
region2 = find( hubspeed > cutin & hubspeed < rated );
region3 = find( hubspeed > rated );

power(region1) =  0;
power(region2) = 0.5 * 1030 * 3.14 * R^2 * hubspeed(region2).^3 * efficiency;
power(region3) = 0.5 * 1030 * 3.14 * R^2 * rated.^3 * efficiency;
AEP = ( Tyear / T ) * trapz( power ) * dt;
AEP_GWh = AEP./1e9;  % 2.2 GWh, within 10% of model
AEP = 2*AEP % two rotors

%% Capacity Factor
ratedpower = mean( power(region3) );
CF = trapz(power) * dt ./ (T * ratedpower) % 0.19


%% speed and power histograms at hubheight

fh9 = figure(9); clf
subplot(1,2,1), hist(hubspeed,[0:.1:4]), xlabel('Speed at hub depth [m/s]'), set(gca,'YLim',[0 1000]), ylabel('10 minute ensembles')
subplot(1,2,2), hist(power./1000,50), xlabel('turbine model output [kW]'),  set(gca,'YLim',[0 1000]), ylabel('10 minute ensembles')
safeExport(fh9, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_hists.png'))

maxspeed = max(hubspeed) % 3.2 m/s

%% correlation with tidal stage 
offset = 89;

fh10 = figure(10); clf
subplot(2,2,1)
plot(t,waterdepth-89,'k'), set(gca,'YLim',[-1.5 3])
set(gca,'fontsize',16,'fontweight','demi'), datetick, ylabel('stage [m]'),

subplot(2,2,3)
plot(t,hubspeed,'k')
set(gca,'fontsize',16,'fontweight','demi'), datetick, ylabel('U_{hub} [m/s]')

subplot(1,2,2)
plot(waterdepth-89, hubspeed,'k.'), axis([-2 4 0 3.5]), grid
set(gca,'fontsize',14,'fontweight','demi'), xlabel('stage [m]'), ylabel('U_{hub} [m/s]')

safeExport(fh10, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_stage.png'))


%% save clean output
fname = fullfile(datadir, 'SeaSpider_Sig250_Rosario.mat');
save(fname, 'u', 'v', 'w', 't', 'z', 'waterdepth', 'hubspeed')

[N bincenters] = hist(hubspeed,[0:.1:4]); output = [bincenters; N;]';
save(fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_HubSpeedHistogram.txt'), 'output', '-ASCII')

%% load raw data to estimate turbulence intensity

DopplerNoise = 0.21 % m/s
plotflag = false;
TurbIntensity = NaN(size(t));

for fi=1:55% 55 total
    load([datadir '/burst/S103607A002_SeaSpider_' num2str(fi) '.mat']) % raw files (1 Hz)
    fh11 = figure(11); clf
    for ti = 1:length(t) 
        thisburst = find( abs( Data.Average_Time - t(ti) ) < median(diff(t))./2  );
        if ~isempty(thisburst) && any(ti ~= ondeck)
            rawhubspeed = ( Data.Average_VelEast( thisburst, hubbin(ti) ).^2 + Data.Average_VelNorth( thisburst, hubbin(ti) ).^2 ).^0.5; 
            rawhubspeed = filloutliers(rawhubspeed,'center');
            TurbIntensity(ti) = ( var(rawhubspeed) - DopplerNoise.^2 ).^0.5 / hubspeed(ti); 
            if plotflag
                figure(11), 
                plot( Data.Average_Time(thisburst), rawhubspeed, 'b.'), hold on
                plot(t(ti),hubspeed(ti), 'ko','linewidth',3,'markersize',14), hold on
            end
        else
        end
    end
    if plotflag
        datetick, ylabel('hub speed [m/s]'), grid
        safeExport(fh11, fullfile(scriptdir, 'TIplots', ['SeaSpider_Sig250_Rosario_file' num2str(fi) '.png']))
    end
end

TurbInensity( imag(TurbIntensity)~= 0 ) = NaN;
TurbInensity = real(TurbIntensity);
TurbInensity( TurbIntensity < 0.01 ) = NaN;

fh12 = figure(12);
binscatter(hubspeed',TurbInensity)
axis([0 4 0 0.5])
xlabel('hub speed [m/s]'), ylabel('Turb. Intensity []'), grid
safeExport(fh12, fullfile(scriptdir, 'SeaSpider_Sig250_Rosario_TI.png'))

hubTI = TurbInensity;
save(fname, 'hubTI', '-append')

function safeExport(fh, filepath)
% Save a figure to filepath as PNG. Prefers exportgraphics, but falls back
% to print() with the painters (vector) renderer if exportgraphics hits
% the known Java rendering bug (HGRasterOutputHelper NullPointerException)
% that can occur on data-dense figures (e.g. large scatter plots). A
% failed exportgraphics call can sometimes leave the figure handle
% invalid, so the fallback is itself guarded and will just skip (with a
% warning) rather than error out the whole script.
    outdir = fileparts(filepath);
    if ~isempty(outdir) && ~isfolder(outdir)
        mkdir(outdir)
    end
    if ~isvalid(fh)
        warning('safeExport:invalidFigure', 'Figure handle is invalid; skipping export of %s', filepath)
        return
    end
    drawnow
    try
        exportgraphics(fh, filepath)
    catch err
        warning('safeExport:exportgraphicsFailed', ...
            'exportgraphics failed (%s); falling back to print() for %s', err.message, filepath)
        if ~isvalid(fh)
            warning('safeExport:figureInvalidAfterFailure', ...
                'Figure became invalid after failed export; skipping %s', filepath)
            return
        end
        try
            print(fh, filepath, '-dpng', '-r150', '-painters')
        catch err2
            warning('safeExport:printFailed', ...
                'print() fallback also failed (%s); skipping %s', err2.message, filepath)
        end
    end
end