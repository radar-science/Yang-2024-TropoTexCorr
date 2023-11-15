%% Test script for the Varying Depth Topo correction for Mogi source.
scrsz = get(0,'ScreenSize');

%% Model setup
X0 = 230;
Y0 = 230;
depth = 730;   % depth below the summit
oX = 0;
oY = 0;
oZ = 0;
aX = 200;
aY = 200;
aZ = 200;
op = 0.30;
nu = 0.25;

%% Prepare X, Y, Z coordinates
wavelength = 0.236;
dem_file = '/Users/yunjunz/Papers/coauthor/2023-Yang-Tropo/data/mintpy/inputs/dem.h5';
inc_angle_file = '/Users/yunjunz/Papers/coauthor/2023-Yang-Tropo/data/mintpy/inputs/incAngle.h5';
dem = double(h5read(dem_file, '/h5'));
dem(dem <750) = 750;
inc_angle = double(h5read(inc_angle_file, '/h5')) * pi / 180;
az_angle = -100.5 * pi / 180;
[len, wid] = size(dem);
[X,Y] = meshgrid(1:wid, 1:len);

X = X - X0;
Y = Y - Y0;

%%
X = reshape(X*20,[1,len*wid]);
Y = reshape(Y*20,[1,len*wid]);
Z = reshape(dem, [1,len*wid]);

%% Predicted displacement without / with considering topo effect
[u,v,w,~] = cdmv(X,Y, Z-depth, oX,oY,oZ,aX,aY,aZ,op,nu);
u = reshape(u, [len,wid]);
v = reshape(v, [len,wid]);
w = reshape(w, [len,wid]);
u(isnan(u)) = 0;
v(isnan(v)) = 0;
w(isnan(w)) = 0;

los = u .* sin(inc_angle) .* sin(az_angle) .* -1 + ...
      v .* sin(inc_angle) .* cos(az_angle) + ...
      w .* cos(inc_angle);

%% Plot
xlims = [min(X), max(X)];
ylims = [min(Y), max(Y)];
% clims = [-wavelength/4, wavelength/4];

figure('Position',[10 scrsz(4)/5 2*scrsz(3)/3 scrsz(4)/3],'Name','Sim Defo at Kirishima');
subplot(2,3,1); imagesc(xlims, ylims, dem);  colorbar;   axis equal; axis tight; title('DEM');
subplot(2,3,2); imagesc(xlims, ylims, u );  colorbar;   axis equal; axis tight; title('East  displacement');
subplot(2,3,3); imagesc(xlims, ylims, v );  colorbar;   axis equal; axis tight; title('North displacement');
subplot(2,3,5); imagesc(xlims, ylims, w );  colorbar;   axis equal; axis tight; title('Up    displacement');
subplot(2,3,6); imagesc(xlims, ylims, w );  colorbar;   axis equal; axis tight; title('LOS   displacement');

colormap jet;
