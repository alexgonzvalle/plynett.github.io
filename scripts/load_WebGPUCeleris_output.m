clf;
clear all;

path_run = 'D:\Development\repos\plynett.github.io\automation\output';

% Load basic parameters
dx = load(fullfile(path_run, 'dx.txt'));
dy = load(fullfile(path_run, 'dy.txt'));
nx = load(fullfile(path_run, 'nx.txt'));
ny = load(fullfile(path_run, 'ny.txt'));

x = [0:1:nx-1] * dx;
y = [0:1:ny-1]' * dy;

% Specify number of time steps to plot
nt = 496;  % need to change manually

% load bathy/topo [m]
fid = fopen(fullfile(path_run, 'bathytopo.bin'), 'r');
bathytopo = fread(fid, [nx, ny], 'float32');
fclose(fid);

time = zeros(nt,1);

for n = 1:nt
    time_fname = [fullfile(path_run, ['time_' num2str(n) '.txt'])];
    time(n) = load(time_fname);  % [sec]

    % free surface elevation [m]
    eta_fname = [fullfile(path_run, ['elev_' num2str(n) '.bin'])]; % Note: Changed from .txt to .bin
    fid = fopen(eta_fname, 'r');
    eta = fread(fid, [nx, ny], 'float32');
    fclose(fid);

%     % x direction flux [m^2/s]
%     P_fname = [fullfile(path_run, ['xflux_' num2str(n) '.bin'])]; % Note: Changed from .txt to .bin
%     fid = fopen(P_fname, 'r');
%     P = fread(fid, [nx, ny], 'float32');
%     fclose(fid);

%     % y direction flux [m^2/s]
%     Q_fname = [fullfile(path_run, ['yflux_' num2str(n) '.bin'])]; % Note: Changed from .txt to .bin
%     fid = fopen(Q_fname, 'r');
%     Q = fread(fid, [nx, ny], 'float32');
%     fclose(fid);

%     % turbulent eddy viscosity from breaking [m/s^2]
%     nu_fname = [fullfile(path_run, ['turb_' num2str(n) '.bin'])]; % Note: Changed from .txt to .bin
%     fid = fopen(nu_fname, 'r');
%     nu = fread(fid, [nx, ny], 'float32');
%     fclose(fid);

    H = eta - bathytopo;  % total water depth
%     u = P ./ H;  % [m/s]
%     v = Q ./ H;  % [m/s]

    drynodes = find(H<=0.0);  % find dry nodes, may need to change 0.0 to some threshold (e.g. 0.1) to eliminate large velocitys in thin flows
    H(drynodes) = NaN;
    eta(drynodes) = NaN;
%     u(drynodes) = NaN;
%     v(drynodes) = NaN;
%     nu(drynodes) = NaN;

%     % to extract properties at particular location
%     xp = 9000;
%     yp = 4000;
%     Hp = interp2(x,y,H',xp,yp);
%     up = interp2(x,y,u',xp,yp);
%     vp = interp2(x,y,v',xp,yp);
%     disp([Hp up vp])

    % Plotting
    figure(1);
    clf
    pcolor(x, y, eta');  % note the transpose '  - the surface data is [nx,ny]
    shading interp;
    axis equal;
    clim([-1 1])
    colorbar;
    title(sprintf('Time = %f seconds', time(n)));
    pause(.01)
end
