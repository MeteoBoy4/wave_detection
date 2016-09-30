clear all;
% for 0.125*0.125 grid, extract grid 20N~40N latitude, 70~100 longitude for
% the main domain
lat = ncread('/run/media/MeteoBoy4/Data/MData/ERA-Interim/2005/div/jan/upmonth/300hpa.nc', 'latitude',401, 161);
lon = ncread('/run/media/MeteoBoy4/Data/MData/ERA-Interim/2005/div/jan/upmonth/300hpa.nc', 'longitude', 561, 241);
div = ncread('/run/media/MeteoBoy4/Data/MData/ERA-Interim/2005/div/jan/upmonth/300hpa.nc', 'd', [561 401 39], [241 161 1]);
vort = ncread('/run/media/MeteoBoy4/Data/MData/ERA-Interim/2005/vorticity/jan/upmonth/300hpa.nc', 'vo', [561 401 39], [241 161 1]);

%The north to south order is annoying, so reorder is needed
lat = lat(end:-1:1);
div = div(:, end:-1:1);
vort = vort(:, end:-1:1);

%Find out the cross-section through specific (longitude, latitude) points,
%and interpolating the div & vort onto that cross section
south = 25; north = 38; west = 81; east = 95;
southlat = find(lat==south, 1, 'first');
westlon = find(lon==west, 1, 'first');
northlat = find(lat==north, 1, 'first');
eastlon  = find(lon==east, 1, 'first');
skew = (northlat - southlat)/(eastlon-westlon);
xs = westlon:1:eastlon;  % the base line use every X coordinate point
ys = skew * (xs - westlon) + southlat; % the line of cross section
divinp = interp2(div, ys, xs);
vortinp = interp2(vort, ys, xs);

% This is the sample's space interval which is used by wcoherence function
% to convert scale to period(wavelength), the interval is considered
% invariant along the specified cross line
duration = sqrt((110/8*(ys(2) - ys(1))) ^ 2 + (110/8 * cos((south + north)/(2*180)*pi)) ^ 2);
% duration = sqrt((110*(ys(2) - ys(1))) ^ 2 + (110 * cos((south + north)/(2*180)*pi)) ^ 2);
% num = size(xs, 2);
% t = 0:duration:duration*(num-1);

% Create the figure of wcoherence, adjust the Y ticklabel to show 'round'
% number, and adjust the three titles. Note: to make use of the
% functionality of converting scales to period, let the duration in space
% represented by seconds
figurew = figure(1);
set(figurew, 'position', [0, 0, 700, 400]);
wcoherence(divinp,vortinp,seconds(duration));
ax = gca;
ax.Title.FontSize = 18;
ytick = round(pow2(ax.YTick), -1);
ax.YTickLabel=ytick;
ylabel(ax, 'Wavelength(km)', 'Fontsize', 12);
xlabel(ax, 'Distance(km)','HorizontalAlignment','center', 'FontSize', 12);
saveas(gcf, 'coherence.pdf');
close(figurew);