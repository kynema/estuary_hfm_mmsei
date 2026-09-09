% Analysis of Rosario Strait ADCP data collected Oct 2024-Jan 2025
% using stablemoor mid-water mooring with up- and down-looking Sig500s
% Velocity and power statistics for OPALCO tidal project
% J. Thomson, Mar 2025

clear all, close all

scriptdir = fileparts(mfilename('fullpath'));
datadir = fullfile(scriptdir, 'data')

load([datadir '/avgd/S100793A017_STBMup_avgd.mat'])
DataUp=Data;
load([datadir '/avgd/S103965A015_STBMdown_avgd.mat'])
DataDown=Data;
clear Data

%% turbine params
cutin = 0.7; % m/s
rated = 2.25; % m/s
efficiency = 0.39;
R = 13.5; % rotor radius (m)
hubdepth = 3.5 + R; % hubheight (m) below surface

%% cell height above seabed (must be dynamic for STBM) for uplooking data

STBMHeight = 45; % m  (45 m at slack tides, when is mooring vertical)
rangebins = Config.Average_BlankingDistance + Config.Average_CellSize/2 + Config.Average_CellSize .* [1:Config.Average_NCells]; % m from instrument (should be same for up- and down-looking)
z = STBMHeight + rangebins; % m above seabed for each cell, when mooring is upright
z = single(z); % change from unit64 to single precision, for later opperations
blowdown = DataUp.Average_Pressure - min(DataUp.Average_Pressure(DataUp.Average_Pressure>40) ); % mooring blow-down in currents
z = ones(length(blowdown),1)*z - blowdown*ones(1,length(z)); % actual height of range bins above seabed for each ensemble

%% rename variables (for convenience)

t = DataUp.Average_Time; % matlab datenum (days from 1-1-0000)
tall = t*ones(1,length(rangebins)); % timestamp as martices to match every bin
u = DataUp.Average_VelEast; % m/s eastward
v = DataUp.Average_VelNorth; % m/s northward
w = ( DataUp.Average_VelUp1 + DataUp.Average_VelUp2) ./ 2; % up
Deltaw = DataUp.Average_VelUp1 - DataUp.Average_VelUp2;
a1 = DataUp.Average_AmpBeam1; a2 = DataUp.Average_AmpBeam2; a3 = DataUp.Average_AmpBeam3; a4 = DataUp.Average_AmpBeam4;
c1 = DataUp.Average_CorBeam1; c2 = DataUp.Average_CorBeam2; c3 = DataUp.Average_CorBeam3; c4 = DataUp.Average_CorBeam4;
STBMdepth = DataUp.Average_Pressure;
waterdepth = 90; % total water depth cannot be determined from STBM pressure, this is surveyed depth


%% remove surface reflections

maxrange = STBMdepth * cosd(20); % geometric max range, very conservative
for i=1:length(t)
    maxbin = find( rangebins > maxrange(i),1);
    if ~isempty(maxbin)
        u(i, maxbin:end) = NaN;
        v(i, maxbin:end) = NaN;
        w(i, maxbin:end) = NaN;
    end
end


%% add in the downlooking data, with consisten extension of z (downwards)

downlookingindices = length(rangebins) + [1:length(rangebins)];
u(:, downlookingindices ) = DataDown.Average_VelEast(1:length(t), :); % m/s eastward
v(:, downlookingindices ) = DataDown.Average_VelNorth(1:length(t), :); % m/s northward
w(:, downlookingindices ) = ( DataDown.Average_VelUp1(1:length(t), :) + DataDown.Average_VelUp2(1:length(t), :)) ./ 2; % up
Deltaw(:, downlookingindices )  = DataDown.Average_VelUp1(1:length(t), :) - DataDown.Average_VelUp2(1:length(t), :);
a1(:, downlookingindices )  = DataDown.Average_AmpBeam1(1:length(t), :);
a2(:, downlookingindices )  = DataDown.Average_AmpBeam2(1:length(t), :);
a3(:, downlookingindices )  = DataDown.Average_AmpBeam3(1:length(t), :);
a4(:, downlookingindices )  = DataDown.Average_AmpBeam4(1:length(t), :);
c1(:, downlookingindices )  = DataDown.Average_CorBeam1(1:length(t), :);
c2(:, downlookingindices )  = DataDown.Average_CorBeam2(1:length(t), :);
c3(:, downlookingindices )  = DataDown.Average_CorBeam3(1:length(t), :);
c4(:, downlookingindices )  = DataDown.Average_CorBeam4(1:length(t), :);


z(:, downlookingindices ) = STBMHeight - ones(length(blowdown),1)*single(rangebins) - blowdown*ones(1,length(rangebins));

seafloor = z < 1;
u(seafloor) = NaN;
v(seafloor) = NaN;
w(seafloor) = NaN;

tall = t*ones(1,length(rangebins)*2); % timestamp as martices to match every bin

% take the merged profiles and make z increase monotonic from seabed
for ti=1:length(t)
    [zsorted zi] = sort( z(ti,:) );
    z(ti,:) = zsorted;
    u(ti,:) = u(ti,zi);
    v(ti,:) = v(ti,zi);
    w(ti,:) = w(ti,zi);
    Deltaw(ti,:) = Deltaw(ti,zi);
    a1(ti,:) = a1(ti,zi);
    a2(ti,:) = a2(ti,zi);
    a3(ti,:) = a3(ti,zi);
    a4(ti,:) = a4(ti,zi);
    c1(ti,:) = c1(ti,zi);
    c2(ti,:) = c2(ti,zi);
    c3(ti,:) = c3(ti,zi);
    c4(ti,:) = c4(ti,zi);
end


%% remove data when out of water

ondeck = find( STBMdepth < 40 );
u(ondeck,:) = NaN;
v(ondeck,:) = NaN;
w(ondeck,:) = NaN;
z(ondeck,:) = NaN;

inwater = find( STBMdepth > 40 ); inwater(1) = []; inwater(end) = [];


%% look at quality metrics / ancillary data

% orientation
fh1 = figure(1); clf
subplot(3,1,1), plot(t(inwater),DataUp.Average_Pitch(inwater), '.', t(inwater),DataUp.Average_Roll(inwater), '.'), datetick, legend('pitch','roll'), set(gca,'YLim',[-30 30])
subplot(3,1,2), plot(t(inwater),DataUp.Average_Heading(inwater), '.'), datetick, legend('heading'), set(gca,'YLim',[0 360])
subplot(3,1,3), plot(t(inwater),DataUp.Average_Temperature(inwater), '.'), datetick, legend('Temperature'), set(gca,'YLim',[8 12])
safeExport(fh1, fullfile(scriptdir, 'STBM_Sig500up_Rosario_PRH.png'))

% amplitude
fh2 = figure(2); clf
for i=1:4
    beamdata = eval(['a' num2str(i) ]);%eval(['DataUp.Average_AmpBeam' num2str(i) '(inwater,:)']);
    subplot(4,1,i), scatter(tall(:), z(:), 2, beamdata(:), 'filled' ),  datetick, hold on, axis tight, set(gca,'YLim',[0 100])
    plot(t, STBMHeight - blowdown + maxrange,'r.', 'linewidth', 1)
    plot(t, STBMHeight - blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
    cb = colorbar; cb.Label.String = 'dB';
end
safeExport(fh2, fullfile(scriptdir, 'STBM_Sig500s_Rosario_Amp.png'))

% correlation
fh3 = figure(3); clf
for i=1:4
    beamdata = eval(['c' num2str(i) ]);
    subplot(4,1,i), scatter(tall(:), z(:), 2, beamdata(:), 'filled' ),  datetick, hold on, axis tight, set(gca,'YLim',[0 100])
    plot(t, STBMHeight - blowdown + maxrange,'r.', 'linewidth', 1)
    plot(t, STBMHeight - blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
    cb = colorbar; cb.Label.String = '%';
end
safeExport(fh3, fullfile(scriptdir, 'STBM_Sig500s_Rosario_Cor.png'))

%  difference in vertical velocity estimates
fh4 = figure(4); clf
scatter(tall(:), z(:), 2, Deltaw(:), 'filled' ), datetick, hold on, axis tight, set(gca,'YLim',[0 100])
plot(t, STBMHeight - blowdown + maxrange,'r.', 'linewidth', 1)
plot(t, STBMHeight - blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
caxis([-0.5 0.5]), cb = colorbar; cb.Label.String = '\Delta w [m/s]';
safeExport(fh4, fullfile(scriptdir, 'STBM_Sig500s_Rosario_Deltaw.png'))


%% component and speed plots (pcolors)

speed = (u.^2 + v.^2 + w.^2).^.5;  % 3D pythagorean thm, in case STBM trim is not perfect.  Does not change AEP result

fh5 = figure(5); clf
colormap(jet)
ax(1)=subplot(3,1,1); scatter(tall(:), z(:), 2, u(:) , 'filled'), datetick, hold on, axis tight, set(gca,'YLim',[0 100])
caxis([-2 2]), cb = colorbar; cb.Label.String = 'east [m/s]';
plot(t, STBMHeight-blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
ax(2)=subplot(3,1,2); scatter(tall(:), z(:), 2, v(:) , 'filled'), datetick, hold on, axis tight, set(gca,'YLim',[0 100])
caxis([-2 2]), cb = colorbar; cb.Label.String = 'north [m/s]';
plot(t, STBMHeight-blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
ax(3)=subplot(3,1,3); scatter(tall(:), z(:), 2, u(:), 'filled' ), datetick, hold on, axis tight, set(gca,'YLim',[0 100])
caxis([-2 2]), cb = colorbar; cb.Label.String = 'up [m/s]';
plot(t, STBMHeight-blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
linkaxes(ax,'x')
safeExport(fh5, fullfile(scriptdir, 'STBM_Sig500s_Rosario_UVW.png'))

fh6 = figure(6); clf
colormap(jet)
scatter(tall(:), z(:), 2, speed(:) , 'filled'), datetick, hold on, axis tight, set(gca,'YLim',[0 100])
caxis([0 3.5]), cb = colorbar; cb.Label.String = 'horizontal speed [m/s]';
plot(t, STBMHeight - blowdown,'k-', 'linewidth', 1), ylabel('z [m]')
safeExport(fh6, fullfile(scriptdir, 'STBM_Sig500s_Rosario_speed.png'))


%% tidal ellipses and u,v scatter
principalaxis = nanmean( atand( v ./ u ) );

fh7 = figure(7); clf
colormap copper
scatter(u(:),v(:), 1, z(:),'filled')
caxis([0 100]), cb = colorbar; cb.Label.String = 'z [m]';
xlabel('u, east [m/s]'), ylabel('v, north [m/s]')
axis([-3 3 -3 3]), grid
safeExport(fh7, fullfile(scriptdir, 'STBM_Sig500s_Rosario_UVscatter.png'))


%% vertical profiles (requires gridding for STBM variable z)
fixedz = [2:81];

speed_fixedz = NaN(length(t),length(fixedz));
power_fixedz =  NaN(length(t),length(fixedz));
axis_fixedz =  NaN(length(t),length(fixedz));

for ti=1:length(t)
    goodpts = isfinite( z(ti,:) + speed(ti,:) );
    if sum(goodpts)>20
        speed_fixedz(ti,:) = interp1( z(ti,goodpts), speed(ti,goodpts), fixedz );
        power_fixedz(ti,:) = interp1( z(ti,goodpts), 0.5*1030*speed(ti,goodpts).^3, fixedz );
        axis_fixedz(ti,:) = interp1( z(ti,goodpts), atand(v(ti,goodpts)./ u(ti,goodpts)), fixedz );
    end
end

fh6b = figure(6); clf
colormap(jet)
pcolor(t, fixedz, speed_fixedz'), shading flat, datetick, hold on, axis tight, set(gca,'YLim',[0 100])
plot(t(inwater), STBMHeight-blowdown(inwater),'k-', 'linewidth', 1), ylabel('z [m]')
caxis([0 3.5]), cb = colorbar; cb.Label.String = 'horizontal speed [m/s]';
safeExport(fh6b, fullfile(scriptdir, 'STBM_Sig500s_Rosario_speed_gridded.png'))


fh8 = figure(8); clf
subplot(1,3,1), plot( nanmean(speed_fixedz), fixedz), hold on, area([0 1.5],[94 94],83), plot([0 1.5],(waterdepth - hubdepth)*[1 1],'k--'),
plot([0 1.5],(waterdepth - hubdepth+R)*[1 1],'k:'), plot([0 1.5],(waterdepth - hubdepth-R)*[1 1],'k:'),
ylabel('z [m]'), xlabel('avg speed [m/s]'),set(gca,'YLim',[0 100])
subplot(1,3,2), plot( nanmean( power_fixedz), fixedz), hold on, area([0 2000],[94 94],83), plot([0 2000],(waterdepth - hubdepth)*[1 1],'k--'),
plot([0 2000],(waterdepth - hubdepth+R)*[1 1],'k:'), plot([0 2000],(waterdepth - hubdepth-R)*[1 1],'k:'),
ylabel('z [m]'), xlabel('avg avail. power desnity [W/m^2]'),set(gca,'YLim',[0 100])
subplot(1,3,3), plot( nanmean(360+axis_fixedz), fixedz), hold on, area([300 320],[94 94],83), plot([300 320],(waterdepth - hubdepth)*[1 1],'k--'),
plot([300 320],(waterdepth - hubdepth+R)*[1 1],'k:'), plot([300 320],(waterdepth - hubdepth-R)*[1 1],'k:'),
ylabel('z [m]'), xlabel('principal axis [deg M]'),set(gca,'YLim',[0 100])
safeExport(fh8, fullfile(scriptdir, 'STBM_Sig500s_Rosario_profiles.png'))



%% AEP estimate for single O2 turbine, based on turbine parameters at top of script,
% note that actual production would be double

% hubheight speed
for i=1:length(t)
    [dz zi] = min( abs( z(i,:) - (waterdepth - hubdepth) ) ) ;
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
AEP_GWh = AEP./1e9;  % 2.45 GWh with 2D speed, 2.47 GWh with 3D speed
AEP = 2*AEP % two rotors

%% Capacity Factor
ratedpower = mean( power(region3) )
CF = trapz(power) * dt ./ (T * ratedpower) % 0.21 with 2D speed, 0.22 with 3D speed


%% speed and power histograms at hubheight

fh9 = figure(9); clf
subplot(1,2,1), hist(hubspeed,[0:.1:4]), xlabel('Speed at hub depth [m/s]'), set(gca,'YLim',[0 1000]), ylabel('10 minute ensembles')
subplot(1,2,2), hist(power./1000,50), xlabel('turbine model output [kW]'),  set(gca,'YLim',[0 1000]), ylabel('10 minute ensembles')
safeExport(fh9, fullfile(scriptdir, 'STBM_Sig500s_Rosario_hists.png'))

fh10 = figure(10); clf
plot(t,hubspeed,'k'), grid
set(gca,'fontsize',16,'fontweight','demi'), datetick, ylabel('U_{hub} [m/s]')
safeExport(fh10, fullfile(scriptdir, 'STBM_Sig500s_Rosario_hubspeed.png'))


maxspeed = max(hubspeed) %

%% save clean output
fname = fullfile(scriptdir, 'STBM_Sig500s_Rosario.mat');
save(fname, 'u', 'v', 'w', 't', 'z', 'hubspeed', '*fixed*')

[N bincenters] = hist(hubspeed,[0:.1:4]); output = [bincenters; N;]';
save(fullfile(scriptdir, 'STBM_Sig500s_Rosario_HubSpeedHistogram.txt'), 'output', '-ASCII')

%% load raw data to estimate turbulence intensity

DopplerNoise = 0.073 % m/s
plotflag = true;
TurbIntensity = NaN(size(t));
flist = dir([datadir '/burst/S100793A017_STBMup_*']);


for fi=1:length(flist)
    load( [datadir '/burst/' flist(fi).name ] ) % raw files (1 Hz)
    fh11 = figure(11); clf
    for ti = 1:length(t)
        thisburst = find( abs( Data.Average_Time - t(ti) ) < median(diff(t))./2  );
        if ~isempty(thisburst) && any(ti ~= ondeck)
            [dz thishubbin ] = min( abs( STBMdepth(ti) - hubdepth - single(rangebins) ) );
            rawhubspeed = ( Data.Average_VelEast( thisburst, thishubbin ).^2 + Data.Average_VelNorth( thisburst, thishubbin ).^2 ).^0.5;
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
        safeExport(fh11, fullfile(scriptdir, 'TIplots', ['STBM_Sig500up_Rosario_file_' flist(fi).name(end-8:end-4) '.png']))
    end
end

TurbInensity( imag(TurbIntensity)~= 0 ) = NaN;
TurbInensity = real(TurbIntensity);
TurbInensity( TurbIntensity < 0.01 ) = NaN;

fh12 = figure(12); clf
binscatter(hubspeed',TurbInensity)
axis([0 4 0 0.5])
xlabel('hub speed [m/s]'), ylabel('Turb. Intensity []'), grid
safeExport(fh12, fullfile(scriptdir, 'STBM_Sig500_Rosario_TI.png'))

hubTI = TurbInensity;
save(fullfile(scriptdir, 'STBM_Sig500s_Rosario.mat'), 'hubTI', '-append')


%% surface waves
% process uplooking Sig500 using sigProcess_all.m from NortekCodes repo
% this calls UVZwaves.m from SWIFTcodes repo
% then saves everything in STBMup_Signature_all_processed.mat

load(fullfile(scriptdir, 'STBM_Sig500up_preprocessed'), 'sigBurst')
Waves = sigBurst;
Waves = rmfield(Waves,'icethickness');
Waves = rmfield(Waves,'icehistogram');
save(fullfile(scriptdir, 'STBM_Sig500s_Rosario.mat'), 'Waves', '-append')


fh21 = figure(21);
subplot(2,1,1), plot([Waves.time], [Waves.sigwaveheight],'.'), set(gca,'YLim',[0 2]), datetick, ylabel('H_s [m]')
subplot(2,1,2), plot([Waves.time], [Waves.peakwaveperiod],'.'), set(gca,'YLim',[0 10]), datetick, ylabel('T_p [s]')
safeExport(fh21, fullfile(scriptdir, 'SBTM_waves_timeseries.png'))

fh22 = figure(22); % wave histogram
binscatter([Waves.peakwaveperiod],[Waves.sigwaveheight])
set(gca,'fontsize',16,'fontweight','demi')
xlabel('Wave T_p [s]')
ylabel('Wave H_s [s]')
grid
safeExport(fh22, fullfile(scriptdir, 'SBTM_waves_hist.png'))

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


% max wave conditions
% [mh mi ] = max([sigBurst.sigwaveheight])
%
% mh =
%
%     1.7899
%
%
% mi =
%
%         4662
%
% sigBurst(mi).peakwaveperiod
%
% ans =
%
%     5.7528
%
%
% L = 9.8*5.7^2 / (6.28)
%
% L =
%
%    50.7010
%
% uorbsurface = 3.14*1.8/5.7
%
% uorbsurface =
%
%     0.9916
%
% uorbhub = uorbsurface * exp(-6.28*hubdepth/L)
%
% uorbhub =
%
%     0.1207