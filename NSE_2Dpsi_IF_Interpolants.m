function  NSE_2Dpsi_IF_Interpolants(Nx,Ny)
% Author: Collin Victor. Last modified on 2020-09-25.
% Adapted from code by Adam Larios
% Solve the 2D Navier-Stokes equations in stream-function form
% Using spectral methods, AB-3 with integrating factor, 2/3 dealiasing.
close all;
% clear all;

% T = 1;
show = 1; %Graphing parameter (plot every $show timesteps)


plot_only_trajectories = false;

% Custom options to save certain things, and to plot/avoid plotting certain
% relevant fields.
%options =  [errorOnly,     saveCheck,  saveData,   saveVars,   makeMovie,  plotSpectrum,   plotObservers,     plotEnstrophy,  u_only]
options =   [0              0           0           1           0           0               0                  0               0];


% da_skip = 50;

warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')


%% =========== Input Parameters ===========

x_min = -pi;
x_max = pi;
y_min = -pi;
y_max = pi;

resolution = 2^10;
switch resolution
    case 2^10
        Nx = 2^10;
        Ny = 2^10;
        T = 10;
        nu = 0.0001; %Viscosity used by Jolly
        alpha = 32/sqrt(Nx^2+Ny^2);
        G = 1e6;%Jolly value
        load('C:\Users\colli\Documents\MATLAB\Research\Desktop\restart files\restart1024.mat','psi_hat_initial')
        dt = 0.005;
    case 2^9
        Nx = 2^9;
        Ny = 2^9;
        T = 1;
        nu = 0.0005;
        alpha = 0.005;
        
        G = 1e5;
        load('.\restart files\restart512_12100.mat','psi_hat_initial')
        dt = 0.01;
        
end





%% === Pre-set options so you can just hit the "run" button in Matlab.
if ~(exist('Nx','var'))
    %     Nx = 2^9;  % Number of x gridpoints (uniform).
    Nx = 2^10;
    %     Nx = 2^7;  % Number of x gridpoints (uniform).
    
elseif mod(Nx,2)
    Nx = Nx+1;
    disp('WARNING: Odd number of gridpoints not supported, adding one to make an even number.');
end
if ~(exist('Ny','var'))
    Ny = Nx; % Number of y gridpoints (uniform).
elseif mod(Nx,2)
    Ny = Ny+1;
    disp('WARNING: Odd number of gridpoints not supported, adding one to make an even number.');
end
if ~(exist('T','var'))
    %     T = 50.0;
    T = 1;
    %     T = 25000;
    %     T = 100;
end
if ~(exist('nu','var'))
    %nu = 0.0001*10^-5;
    %     nu = .9;
    %     nu = 0.00001;
    nu = 0.0001; %Viscosity used by Jolly
    % nu = 0;
    %     nu = 0.005;
    
    
end
if ~(exist('alpha','var'))
    alpha = 32/sqrt(Nx^2+Ny^2);
    %     alpha = 0;
    %     alpha = 1;
end
if ~(exist('G','var'))
    %     G = 2.0; % Grashof number
    %     G = 5e5;
    %         G = 2.5e6; % Gesho-Olson-Titi value
    G = 1e6;%Jolly value
    %    15.2048257930675 if we use little l1 norm, or 10.7546058538537
    
end
if ~(exist('U','var'))
    %     U = 1.0; % Velocity L^2 norm
end

% % Seed RNG
% % RNG_seed = randi(10000);
% RNG_seed = 12345;
% randn('state',RNG_seed);

%% =========== Grid Set Up ===========
Lx = x_max - x_min;
Ly = y_max - y_min;
L = sqrt(Lx*Ly);

p.Lx = Lx;
p.Ly = Ly;
p.L = L;

p.parseval = sqrt(Lx*Ly)/Nx/Ny; % Multiply by norm(u_hat(:)) to get L^2 norm

% Physical space
dx = (x_max - x_min)/Nx; % Assume uniform grid.
x  = x_min + (0:(Nx-1))*dx;
dy = (y_max - y_min)/Ny; % Assume uniform grid.
y  = y_min + (0:(Ny-1))*dy;

p.dx = dx;
p.dy = dy;
p.x = x;
p.y = y;
[p.X,p.Y] = meshgrid(x,y);


x_pts = reshape(p.X, Nx^2,1);
y_pts = reshape(p.Y, Ny^2,1);
p.coordinates = [x_pts, y_pts];
% x_pts_per = [];
% y_pts_per = [];
% for(i = -1:1)
%     for(j=-1:1)
%         x_pts_per = [x_pts_per; x_pts + Lx*i];
%         y_pts_per = [y_pts_per; y_pts + p.Ly*j];
%     end
% end

% p.coordinates_periodic = [x_pts_per,y_pts_per];



% Wave numbers k (we also multiply k by i for convenience).
kx    =      ([0:Nx/2, -Nx/2+1:-1]*(2*pi/Lx)).';
p.ikx = (1i*[0:((Nx/2)-1) 0 ((-Nx/2)+1):-1]*(2*pi/Lx)).';
p.kx = ([0:((Nx/2)-1) 0 ((-Nx/2)+1):-1]*(2*pi/Lx)).';
kx_sq =    ([0:Nx/2, -Nx/2+1:-1].^2*(2*pi/Lx)^2).';

% Wave numbers k (we also multiply k by i for convenience).
ky    =    [0:Ny/2, -Ny/2+1:-1]*(2*pi/Ly);
p.iky = 1i*[0:((Ny/2)-1) 0 ((-Ny/2)+1):-1]*(2*pi/Ly);
p.ky = [0:((Ny/2)-1) 0 ((-Ny/2)+1):-1]*(2*pi/Ly);
ky_sq =    [0:Ny/2, -Ny/2+1:-1].^2*(2*pi/Ly)^2;

% Fourier Truncation of high wave modes
% p.trunc_num = 50; %Arbitrary choice of truncation numer
% p.trunc_array = ones(Nx,Ny);
% for j = 1:Ny
%     for i = 1:Nx
%         if(p.kx(i)^2 + p.ky(j)^2 >= p.trunc_num^2)
%             p.trunc_array(i,j) = 0;
%         end
%     end
% end

% spy(p.trunc_array)

% Wave numbers of Laplacian operator.
p.k_lap     = zeros(Nx,Ny);
p.k_lap_inv = zeros(Nx,Ny);
p.k1sq_m_k2sq = zeros(Nx,Ny);
p.helm_inv = zeros(Nx,Ny);
for j = 1:Ny
    for i = 1:Nx
        norm_k_sq = (kx_sq(i) + ky_sq(j));
        p.k_lap(i,j) = -norm_k_sq;
        p.k_lap_inv(i,j) = -1/norm_k_sq;
        
        p.k1sq_m_k2sq(i,j) = kx_sq(i)-ky_sq(j);
        p.k1k2(i,j) = kx(i)*ky(j);
        
        p.helm_inv(i,j) = 1./(1+alpha^2*norm_k_sq);
    end
end
p.k_lap_inv(abs(p.k_lap_inv) == Inf) = 0;
p.k_lap_inv(isnan(p.k_lap_inv)) = 0;
p.k_lap_inv(1,1) = 0;

if (alpha<sqrt(eps))
    p.helm_inv = 1;
end

p.Nx    = Nx;
p.Ny    = Ny;
p.Nx_sq = Nx^2;
p.Ny_sq = Ny^2;

% Make a mask for dealiasing
% dealias_modes_x = ceil(Nx/3):(Nx - floor(Nx/3) + 1);
% dealias_modes_y = ceil(Ny/3):(Ny - floor(Ny/3) + 1);
dealias_mask = zeros(Nx,Ny);
dealias_mask(ceil(Nx/3):(Nx - floor(Nx/3) + 1),:) = 1;
dealias_mask(:,ceil(Ny/3):(Ny - floor(Ny/3) + 1)) = 1;
p.dealias_modes = [1;find(dealias_mask == 1)]; % Also enforce mean-zero
clear dealias_mask;

% p.alpha = alpha;
p.nu = nu;
p.ti = 0;


%% =========== Initial Conditions ===========
% Commment/Uncomment the initial condition you want or leave alone to use
% restart file loaded above.

% % --- psi Bump Function ---
% psi = zeros(Nx,Ny);
% for j = 1:Ny
%     for i = 1:Nx
%         psi(i,j) = (exp(1/(x(i)^2 + y(j)^2 - 1))).*(x(i)^2 + y(j)^2 < 1);% + 0.01*rand(1);
%     end
% end
% psi(abs(psi) == Inf) = 0; % 0*Inf = Inf, so zero-out those if they arise.
% psi(isnan(psi)) = 0; % 0*NaN = NaN, so zero-out those too.
% % ---

% % --- psi Taylor-Green ---
% psi = zeros(Nx,Ny);
% for j = 1:Ny
%     siny = sin(y(j));
%     for i = 1:Nx
%         psi(i,j) = sin(x(i)).*siny;
%     end
% end
% % ---

% % --- psi Annular region with random coefficients ---
% psi_hat = zeros(Nx,Ny);
% inner_wave_number = 2;
% outer_wave_number = 4;
% inner_radius_sq = inner_wave_number^2*(2*pi/Lx)^2;
% outer_radius_sq = outer_wave_number^2*(2*pi/Lx)^2;
%
% for j = 1:(Ny/2)
%     for i = [1:(Nx/2), (Nx/2 + 2):Nx] % Jump over Nyquist frequency
%          if (abs(p.k_lap(i,j)) < outer_radius_sq) && ...
%                  (abs(p.k_lap(i,j)) > inner_radius_sq)
%             psi_hat(i,j) = randn(1) + 1i*randn(1);
%          end
%     end
% end
% psi_hat = fft2(ifftn(psi_hat,'symmetric'));
% psi_hat(abs(psi_hat) < sqrt(eps)) = 0; % Kill off errors
%
% % Normalize in L^2 norm
% psi_hat = psi_hat*(U/(norm(psi_hat)*p.parseval));
%
% psi = ifftn(psi_hat);
%
% norm(ifftn(psi_hat))*sqrt(dx*dy)




% psi_hat_initial = zeros(Nx,Ny);
psi_hat = psi_hat_initial.*1;
psi = ifftn(psi_hat);

% ---

%% =========== Constant Forcing ===========

p.f_hat = OlsonTitiForcing(Nx,Ny,nu,G);

%% =========== Initialization and Setup ===========
psi_hat = fftn(psi);


omega_hat = -p.k_lap.*psi_hat;
omega = ifftn(omega_hat,'symmetric');

% Compute velocity norm for CFL and for initialization
en = 0;
for j = 1:p.Ny
    en = en + norm(p.iky(j)*psi_hat(:,j))^2;
    for i = 1:p.Nx
        en = en + norm(p.ikx(i)*psi_hat(i,:))^2;
    end
end
en = sqrt(en)/Nx/Ny; %#ok<*NASGU>

p.dt = dt;

t = 0:dt:T;
nT_steps = length(t);
p.nT_steps  = nT_steps;

p.t = t;
p.T = T;


energy    = zeros(1,nT_steps);
enstrophy = ones(1,nT_steps)*NaN;
enstrophy_Linf = ones(1,nT_steps)*NaN;
L2H1_sq   = zeros(1,nT_steps);


E  = exp(nu*p.k_lap.*dt/2);
E2 = exp(p.nu*p.k_lap.*dt);
E3 = exp(p.nu*p.k_lap.*2*dt);
E4 = exp(p.nu*p.k_lap.*3*dt);

past_rhs1 = [];
past_rhs2 = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Assimilation %%
mu = 500;
p.mu = mu;

loadvars = DataAssimilationVariables_NSE(p);
% loadvars = [];
p.size_vars = length(loadvars);
% If there are no v's to run, it should only run u.
if(isfield(loadvars,'v')==0)
    vars = [];
    p.size_vars = 0;
end

vars = initialize_default_var(loadvars,p);


plotvar = vars;


%% Plots
if(plot_only_trajectories)
    varfigures = [];
    varplots = [];
    
    baseplots = [];
    basefigures = [];
    basefigures = figure;
    baseplots(1) = pcolor(x,y,omega);
    axis('square');
    axis tight;
    colormap jet; % Bone, copper, jet
    shading interp; % flat, interp, faceted
    lighting phong;
    colorbar;
    title(sprintf('Reference Vorticity at t = %1.2f',0));
    drawnow;
    
    hold on;
    for i = 1:p.size_vars
        scatter(vars(i).nodes_coordinates(:,2),vars(i).nodes_coordinates(:,1),100,'x', 'black', 'LineWidth', 2);
    end
    drawnow;
    
end


if(~plot_only_trajectories && ~options(1))
    varfigures = [];
    varplots = [];
    
    if(p.size_vars == 0) %#ok<*ISMT>
        baseplots = [];
        basefigures = [];
        basefigures = figure;
        baseplots(1) = pcolor(x,y,omega);
        axis('square');
        axis tight;
        colormap jet; % Bone, copper, jet
        shading interp; % flat, interp, faceted
        lighting phong;
        colorbar;
        %        caxis([-.25,.25]);
        title(sprintf('Reference Vorticity at t = %1.2f',0));
        drawnow;
        
        
        
        
    end
    
    [u1,u2] = psi_converter(psi_hat, p);
    
    for i=1:p.size_vars
        %Initialize interpolated v for plotting.
        
        
        if(vars(i).interp_type~="HOT")
            %Produce interpolated u observations
            u_obs = compute_observations(u1, u2, vars(i), p);
            
            %Produce v interpolation
            v_interp = compute_Ih_v(vars(i),p);
            
            vars(i).v1_mu = u_obs - v_interp;
        else
            vars(i).v1_mu = psi_hat.*vars(i).trunc_array;
            
        end
        
        figname = sprintf('Observer Type: %s, Smoothing Type: %s', vars(i).observer_type,  vars(i).interp_type);
        varfigures(i) = figure('Position',[1 1 1100 500],'Color',[1 1 1], 'Name', figname); %#ok<*AGROW>
        %                 subplot(2,2,1);
        %                         varfigures(i) =  tightfig(varfigures(i));
        %
        subplot(2,2,1);
        varplots(i,1) = pcolor(x,y,omega);
        axis('square');
        axis tight;
        colormap jet; % Bone, copper, jet
        shading interp; % flat, interp, faceted
        lighting phong;
        colorbar;
        %         caxis([-2,2]);
        axmin = max(max(max(omega)),1);
        %         caxis([-axmin axmin]);
        
        %         caxis([-.25,.25]);
        title(sprintf('Reference Vorticity at t = %1.2f',0));
        
        %                 title(sprintf('Vorticity at t = %1.2f',t(1)));
        %     get(p1,'cdata');
        set(varplots(i,1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
        %         drawnow;
        hold on;
        if(options(7))
            varplots(i,5) = scatter(vars(i).nodes_coordinates(:,2),vars(i).nodes_coordinates(:,1), 'filled');
        end
        
        
        subplot(2,2,2);
        %         intdiff = DA_interp(psi_hat, vars(i), p);
        
        %         intdiff = intdiff/(norm(intdiff(:))*sqrt(Lx*Ly)/Nx/Ny);
        varplots(i,2) = pcolor(x,y,omega);
        set(varplots(i,2),'cdata',ifftn( -p.k_lap.*vars(i).v1_mu,'symmetric'));
        axis('square');
        %         colormap jet;
        colormap jet; % Bone, copper, jet
        shading interp; % flat, interp, faceted
        lighting phong;
        colorbar;
        %         caxis([-10,10]);
        
        title(sprintf('Interpolated Vorticity at t = %1.2f',0));
        
        %         hold on;
        %         varplots(i,5) = scatter(vars(i).nodes_coordinates(:,2),vars(i).nodes_coordinates(:,1),5,'filled');
        %
        
        
        %         vars(i).i_nodes_x;
        %         scatter3(X,Y,Z);
        %         stem3(X,Y,V,'.','color','k','MarkerSize',15)
        
        
        subplot(2,2,3);
        %         varplots(i,3) = semilogy(0:dt:T,vars(i).error,'--','LineWidth',1);
        varplots(i,3) = semilogy(0:dt:T,vars(i).ens_umv,'--','LineWidth',1);
        title('Error for Simulated Solution');
        xlabel('Time');
        ylabel('Enstrophy of u-v');
        %         axis([0 T 1e-20 1e-0]);
        %         hold on;
        %         text = sprintf('vars(%d).error',e);
        %         set(errorPlots(e), 'YDataSource',text);
        %         varplots(i,3) = pcolor(x,y,vars(i).v - ifftn( -p.k_lap.*psi_hat,'symmetric'));
        %         hold on;
        %         axis('square');
        %         % colormap copper; % Bone, copper, jet
        %         shading interp; % flat, interp, faceted
        %         lighting phong;
        %         %         caxis([-3,3]);
        % %         title(sprintf('Vorticity Differences at t = %1.2f',t(1)));
        %         get(varplots(i,3),'cdata');
        %         caxis([0,.25]);
        %         title(sprintf('Absolute Vorticity Difference at t = %1.2f',0));
        
        
        %         drawnow;
        
        subplot(2,2,4);
        varplots(i,4) = pcolor(x,y,vars(i).v);
        hold on;
        axis('square');
        % colormap copper; % Bone, copper, jet
        shading interp; % flat, interp, faceted
        lighting phong;
        %         caxis([-3,3]);
        title(sprintf('Vorticity at t = %1.2f',t(1)));
        get(varplots(i,4),'cdata');
        colorbar;
        %         caxis([-.25,.25]);
        title(sprintf('Simulated Vorticity at t = %1.2f',0));
        
        drawnow;
        
        
    end
    if(options(5))
        
        movie1_filename = 'ComparisonMovie.gif';
        
        
        % Capture the plot as an image
        frame = getframe(varfigures(1));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Write to the GIF File
        imwrite(imind,cm,movie1_filename,'gif','DelayTime',0.1,  'Loopcount',inf);
        
    end
    if(options(6)) %Plot Spectrum
        spectrum = zeros(1,Nx/2);
        %         hold off;
        for j = 1:(Ny/2)
            for i = 1:(Nx/2)
                wave_index = floor(sqrt(i^2 + j^2));
                if wave_index <= Nx/2
                    spectrum(wave_index) = spectrum(wave_index) + abs(psi_hat(i,j))^2;
                end
            end
        end
        spectralPlot = figure;
        spectralPlotData = loglog((1:Nx/2),sqrt(spectrum)/Nx^2,'b','LineWidth',1);
        hold on;
        loglog((1:Nx/2),1e-2*(1:Nx/2).^(-5/3),'k','LineWidth',1);
        loglog((1:Nx/2),1e2*(1:Nx/2).^(-5),'color',[0 .5 0],'LineWidth',1);
        loglog([Nx/3,Nx/3],[1e-30,1e2],'r','LineWidth',1);
        title(sprintf('Spectrum at t = %1.2f',t(1)));
        
        for vi = 1: p.size_vars
            spectrum_var = zeros(1,Nx/2);
            for j = 1:(Ny/2)
                for i = 1:(Nx/2)
                    wave_index = floor(sqrt(i^2 + j^2));
                    if wave_index <= Nx/2
                        %                         spectrum_var(wave_index) = spectrum_var(wave_index) + abs(psi_hat(i,j) - vars(vi).v_hat(i,j))^2;
                        spectrum_var(wave_index) = spectrum_var(wave_index) + abs(vars(vi).v_hat(i,j))^2;
                        
                    end
                end
            end
            varplots(vi,5)= loglog((1:Nx/2),sqrt(spectrum_var)/Nx^2,'-');
        end
        
        axis('tight');
        legendInfo = cell(1,p.size_vars + 4);
        legendInfo{1} = 'Spectrum';
        [legendInfo{2:4}] = deal('r^{-5/3}','r^{-5}','Resolution Cutoff');
        for(vi=1:p.size_vars)
            legEntry = sprintf("Obs: %s, Interp: %s",vars(vi).observer_type, vars(vi).interp_type);
            %             switch vars(i).interp_type
            %                 case "Standard"
            %                     legEntry = 'Standard';
            %                 case "Thin Sweep"
            %                     legEntry = 'Sweeps';
            %                 case "Thin Sweep (x2)"
            %                     legEntry = 'Sweeps';
            %                 case "Random Walkers"
            %                     legEntry = 'Creeps';
            %                 case "Randoms"
            %                     legEntry = 'Bleeps';
            %                 case "Thick Sweeps"
            %                     legEntry = 'Sweeps (var 1)';
            %                 case "Airplane"
            %                     legEntry = 'Airplane';
            %                 case "Lagrangian"
            %                     legEntry = 'Lagrangian'
            %             end
            legendInfo{vi+4} = legEntry;
        end
        legend(legendInfo);
        drawnow;
    end
    if(options(8))
        ensFig = figure;
        ensPlot = semilogy(0:dt:T,enstrophy,'--','LineWidth',1);
        title('Enstrophy Plot');
        xlabel('Time');
        ylabel('Enstrophy');
        
    end
    
end

if(options(1))
    wb = waitbar(0,'Running the main loop.');
end


%end_flag ends program end early if convergence or blowup is detected in
%all v's
if(~isempty(vars))
    end_flag = zeros(1,p.size_vars);
end



%% Main Loop %%
for ti = 1:nT_steps
    
    p.ti = p.t(ti);
    %% Plotting
    if (plot_only_trajectories && mod(ti-1,show)==0)
        for i = 1:p.size_vars
            scatter(vars(i).nodes_coordinates(:,2),vars(i).nodes_coordinates(:,1),1,'black', 'filled', 'LineWidth',1.5 );
        end
        title(sprintf('Reference Vorticity at t = %1.2f',t(ti)));
        drawnow;
    end
    if(options(1))
        waitbar(ti/nT_steps,wb,'Running the main loop.');
    end
    if (~plot_only_trajectories && ~options(1)&&mod(ti-1,show)==0)
        if(p.size_vars==0)
            figure(basefigures);
            set(baseplots(1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
            title(sprintf('Reference Vorticity at t = %1.2f',t(ti)));
            colorbar;
            drawnow;
            
            
        end
        
        for(i = 1:p.size_vars)
            figure(varfigures(i));
            subplot(2,2,1);
            set(varplots(i,1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
            title(sprintf('Reference Vorticity at t = %1.2f',t(ti)));
            colorbar;
            if(options(7))
                set(varplots(i,5),'ydata',vars(i).nodes_coordinates(:,1));
                set(varplots(i,5),'xdata',vars(i).nodes_coordinates(:,2));
            end
            
            subplot(2,2,2);
            axmax = max(max(abs(ifftn( -p.k_lap.*psi_hat,'symmetric'))));
            set(varplots(i,2),'cdata',ifftn( -p.k_lap.*vars(i).v1_mu,'symmetric'));
            title(sprintf('Interpolated Vorticity at t = %1.2f',t(ti)));
            colorbar;
            
            
            subplot(2,2,3);
            varplots(i,3) = semilogy(0:dt:T,vars(i).ens_umv,'-','LineWidth',1);
            hold on;
            semilogy(0:dt:T,vars(i).error,'-','LineWidth',1);
            hold off;
            legend('Enstrophy Error', 'Psi L^2 Error');
            
            title('Error of Simulated Solution');
            xlabel('Time');
            ylabel('Error');
            
            subplot(2,2,4);
            set(varplots(i,4),'cdata',( ifftn(-p.k_lap.*fftn(vars(i).v), 'symmetric' )));
            title(sprintf('Simulated Vorticity at t = %1.2f',t(ti)));
            colorbar;
            
            
            drawnow;
            
        end
        if(options(5))
            %          Capture the plot as an image
            frame = getframe(varfigures(1));
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            % Write to the GIF File
            imwrite(imind,cm,movie1_filename,'gif','DelayTime',0.1, 'WriteMode','append');
        end
        if(options(6))
            %             hold off;
            spectrum = zeros(1,Nx/2);
            for j = 1:(Ny/2)
                for i = 1:(Nx/2)
                    wave_index = floor(sqrt(i^2 + j^2));
                    if wave_index <= Nx/2
                        spectrum(wave_index) = spectrum(wave_index) + abs(psi_hat(i,j))^2;
                    end
                end
            end
            figure(spectralPlot);
            set(spectralPlotData,'ydata',sqrt(spectrum)/Nx^2);
            
            
            for vi = 1: p.size_vars
                spectrum_var = zeros(1,Nx/2);
                for j = 1:(Ny/2)
                    for i = 1:(Nx/2)
                        wave_index = floor(sqrt(i^2 + j^2));
                        if wave_index <= Nx/2
                            spectrum_var(wave_index) = spectrum_var(wave_index) + abs(vars(vi).v_hat(i,j))^2;
                            
                        end
                    end
                end
                set(varplots(vi,5), 'ydata', sqrt(spectrum_var)/Nx^2);
            end
            
            
            title(sprintf('Spectrum at t = %1.2f',t(ti)));
            drawnow;
        end
        if(options(8))
            figure(ensFig);
            ensPlot = semilogy(0:dt:T,enstrophy,'-','LineWidth',1);
            
        end
    end
    
    
    %% Update with Data Assimilation
    
    if(~plot_only_trajectories)
        [u1,u2] = psi_converter(psi_hat, p);
        for(i = 1: p.size_vars)
            if(end_flag(i)==0)
                if(vars(i).observer_type == "Lagrangian")
                    if(t>vars(i).resetTime)
                        vars(i) = reset(vars(i));
                    end
                end
                v_hat = fft2(vars(i).v);
                t_start = cputime;
                
                [~,v1] = compute_rhs_hat(t(ti)          ,v_hat     ,p);
                
                vars(i).v1_mu = zeros(Nx,Ny);
                %             if(vars(i).interp_type == "HOT")
                %Produce interpolated u observations
                u_obs = compute_observations(u1, u2, vars(i), p);
                
                %Produce v interpolation
                v_interp = compute_Ih_v(vars(i),p);
                
                vars(i).v1_mu = u_obs - v_interp;
                
                v1 = v1 +  vars(i).mu*vars(i).v1_mu; %Data is interpolated at stream function level
                v2 = vars(i).v2;
                v3 = vars(i).v3;
                %         if(isempty(vars(i).v2)||isempty(vars(i).v3))
                %             v_hat = E2.*v_hat + dt*E2.*v1;
                %         else
                
                
                
                
                %                 v_hat = E2.*v_hat + dt*E2.*v1;
                %         v_hat = E2.*(v_hat + dt*v1);
                % AB3 with IF
                % v_(n+1) = exp(dt*L)v_n
                %           + dt(23/12) exp(dt*L)  N(v_n,t_n)
                %           + dt(-4/3)  exp(2dt*L) N(v_(n-1),t_(n-1))
                %           + dt(5/12)  exp(3dt*L) N(v_(n-2),t_(n-2))
                
                v_hat = E2.*v_hat + dt*(23/12)*E2.*v1 + dt*(-4/3)*E3.*v2 + dt*(5/12)*E4.*v3;
                
                %Euler Time stepping with IF
                %             v_hat = E2.*v_hat + dt*E2.*v1;
                %         end
                vars(i).v3 = vars(i).v2;
                vars(i).v2 = v1;
                %         vars(i).v2 = v_hat;
                vars(i).v_hat = v_hat;
                vars(i).v = ifftn(v_hat, 'symmetric');
                vars(i) = updateNodes(vars(i),p,psi_hat);
            end
            t_end = cputime;
            vars(i).cpu_time(ti+1) = vars(i).cpu_time(ti) + t_end - t_start;
            %         plotvar(i) = updateNodes(plotvar(i),p);
            plotvar(i) = vars(i);
        end
    end
    
    
    [u1,u2] = psi_converter(psi_hat, p);
    for(i = 1: p.size_vars)
        vars(i) = updateNodes(vars(i),p,psi_hat);
    end
    
    %% Update u in main loop %%
    if(isempty(past_rhs1) || isempty(past_rhs2))
        
        %         Update with RK4 for integrating Factors:
        [en,k1] = compute_rhs_hat(t(ti)         ,psi_hat                 ,p); %#ok<*ASGLU>
        [~ ,k2] = compute_rhs_hat(t(ti) + 0.5*dt,E.*(psi_hat + 0.5*dt*k1),p);
        [~ ,k3] = compute_rhs_hat(t(ti) + 0.5*dt, E.*psi_hat + 0.5*dt*k2 ,p);
        [~ ,k4] = compute_rhs_hat(t(ti) +     dt,E2.*psi_hat +  dt*E.*k3 ,p);
        psi_hat = E2.*psi_hat + (dt/6)*(E2.*k1 + 2*E.*(k2 + k3) + k4);
        
        
        past_rhs2 = past_rhs1;
        past_rhs1 = k1;
        
        psi = ifftn(psi_hat, 'symmetric');
        
    else
        
        % Update with AB3 for integrating Factors:
        [en,k1] = compute_rhs_hat(t(ti)             , psi_hat                ,p);
        psi_hat = E2.*(psi_hat + dt*(23/12*k1)) + E3.*(-4/3*dt*past_rhs1) + E4.*(5/12*dt*past_rhs2);
        psi = ifftn(psi_hat, 'symmetric');
        
        past_rhs2 = past_rhs1;
        past_rhs1 = k1;
    end
    
    
    
    omega_hat_u = p.k_lap.*psi_hat;
    omega_u = real(ifftn(omega_hat_u));
    enstrophy(ti) = norm(omega_hat_u)*p.parseval;
    enstrophy_Linf(ti) = max(max(omega_u));
    
    %% Calculate Error %%
    for(i = 1:p.size_vars)
        errorcheck = 1;
        if(end_flag(i)==0)
            errorcheck = (sqrt(dx*dy)*norm((vars(i).v(:) - psi(:))));
            omega_hat_v = p.k_lap.*vars(i).v_hat;
            omega_v = real(ifftn(omega_hat_v));
            vars(i).ens_u(ti) = enstrophy(ti);
            vars(i).ens_u_Linf(ti) = enstrophy_Linf(ti);
            vars(i).ens_v(ti) = norm(omega_hat_v)/Nx/Ny;
            vars(i).ens_v_Linf(ti) = max(max(omega_v));
            vars(i).ens_umv(ti) = norm(omega_hat_u - omega_hat_v)*p.parseval;
            vars(i).ens_umv_Linf(ti) = max(max(omega_u - omega_v));
            if(errorcheck > 1 || errorcheck < 1e-14)
                end_flag(i) = 1;
            end
            vars(i).error(ti) = errorcheck;
        end
    end
    if(~p.size_vars == 0)
        if(sum(end_flag) == p.size_vars)
            break;
        end
    end
    
    
    
end

%% MISC. %%
if(options(1))
    delete(wb);
end

for(i=1:p.size_vars)
    vars(i).cpu_time(end) = [];
end


if(plot_only_trajectories)
    hold on;
    set(baseplots(1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
    axis('square');
    axis tight;
    colormap jet; % Bone, copper, jet
    shading interp; % flat, interp, faceted
    lighting phong;
    colorbar;
    %        caxis([-.25,.25]);
    title(sprintf('Reference Vorticity at t = %1.2f',T));
    drawnow;
    
    for i = 1:p.size_vars
        scatter(vars(i).nodes_coordinates(:,2),vars(i).nodes_coordinates(:,1),100,'black','filled','s');
    end
    
    
    set(gca, 'FontSize',24);
    drawnow;
    %     f = baseplots(1);
    %     aH = axes;
    %     aHPos = aH.Position;
    %     widthDummy = 0.1;
    %     dummyPos = [aHPos(1)+aHPos(3), aHPos(2), widthDummy, aHPos(4)];
    %     dummyAH = axes('Position', dummyPos,'Color',f.Color, 'XColor',f.Color,...
    %                 'YColor',  f.Color,'XTick',[], 'YTick', []);
end


%% Plotting Error %%
figure;
legendInfo = cell(1,p.size_vars);
for(i=1:p.size_vars)
    semilogy(0:dt:T,vars(i).error,'-','LineWidth',1);
    hold on;
    legEntry = 'Other';
    switch vars(i).observer_type
        case "Uniform"
            legEntry = 'Uniform Grid';
        case "Thin Sweep"
            legEntry = 'Sweeps';
        case "Thin Sweep (x2)"
            legEntry = 'Sweeps';
        case "Random Walkers"
            legEntry = 'Creeps';
        case "Randoms"
            legEntry = 'Bleeps';
        case "Thick Sweeps"
            legEntry = 'Sweeps (var 1)';
        case "Airplane"
            legEntry = 'Airplane';
        case "Lagrangian"
            legEntry = 'Lagrangian';
        case "Local"
            legEntry = 'Local';
    end
    legendInfo{i} = legEntry;
end
xlabel('Time (simulated)');
ylabel('L^2 norm of u-v');
title('Error for Various Methods');
legend(legendInfo);



figure;
% legendInfo = cell(1,p.size_vars);
for(i=1:p.size_vars)
    n_time=numel(vars(i).cpu_time);
    
    semilogy(vars(i).cpu_time,vars(i).error(1:n_time),'-','LineWidth',1);
    hold on;
end
xlabel('Time (cpu)');
ylabel('L^2 norm of u-v');
title('Error for Various Methods');
legend(legendInfo);






if (options(3) == 1)
    save('psi_hat.mat','psi_hat');
end

if(options(4))
    %     date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
    for(i = 1: p.size_vars)
        date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
        var_type = 'Misc.';
        switch vars(i).observer_type
            case "Uniform"
                var_type = 'UniformGrid';
            case "Thin Sweeps"
                var_type = 'Sweeps';
            case "Random Walkers"
                var_type = 'Creeps';
            case "Randoms"
                var_type = 'Bleeps';
            case "Lagrangian"
                var_type = 'Lagrange';
            case "Thick Sweeps"
                var_type = 'SweepsOnBox';
                
        end
        var = vars(i);
        save(sprintf('Data/%s_vars_%s.mat',date_string,var_type),'var', 'p');
        pause(1);
    end
    %    save('vars.mat','vars');
end

% subplot(1,2,2);
% plot(t,enstrophy);
% title('enstrophy');
% ylim([0,max(enstrophy)]);
%
% % Make a sound when done.
% sound(2*sin(2*pi*300*(0:1/20500:0.1)))

end % end main function

function [var] =  updateNodes(var, p, psi_hat)
switch var.observer_type
    case {"Thin Sweeps" "Thin Sweeps (x2)"}
        %     var.i_nodes_y
        new_indices_y = floor(var.i_nodes_y + var.velocity_y*p.dt/p.dy);
        new_indices_x = floor(var.i_nodes_x + var.velocity_x*p.dt/p.dx);
        
        var.i_nodes_y = var.indexing_array(new_indices_y);
        var.i_nodes_x = var.indexing_array(new_indices_x);
        
        %     var.i_nodes_y = sort(mod(var.i_nodes_y+var.velocity, p.Ny+1));
        %     if(var.i_nodes_y(1) == 0)
        %         [~,index] = max(var.i_nodes_y(2:end)-var.i_nodes_y(1:end-1));
        %         if(maximum > 1)
        %          var.i_nodes_y(1)
        %             var.i_nodes_y(1:index) = var.i_nodes_y(1:index)+1;
        %             var.i_nodes_y = (v.x_nodes+p.dx)./p.dx;
        %     end
    case {"Thick Sweeps"}
        new_indices_x = floor(var.i_nodes_coordinates(:,1) + var.velocity_x*p.dt/p.dx);
        new_indices_y = floor(var.i_nodes_coordinates(:,2) + var.velocity_y*p.dt/p.dy);
        var.i_nodes_coordinates(:,1) = var.indexing_array(new_indices_x);
        var.i_nodes_coordinates(:,2) = var.indexing_array(new_indices_y);
        
        
        new_zeros_x = floor(var.zeros_coordinates(:,1) + var.velocity_x*p.dt/p.dx);
        new_zeros_y = floor(var.zeros_coordinates(:,2) + var.velocity_y*p.dt/p.dy);
        var.zeros_coordinates(:,1) = var.indexing_array(new_zeros_x);
        var.zeros_coordinates(:,2) = var.indexing_array(new_zeros_y);
        
        var.nodes_coordinates = [p.x(var.i_nodes_coordinates(:,1))', p.y(var.i_nodes_coordinates(:,2))'];
        
        
    case "Random Walkers"
        for i = 1:length(var.i_nodes_coordinates)
            %         return
            switch randi(5) % Pick a random number form {1,2,3,4}
                case 1 % Step North
                    v1 = 0; v2 = 1;
                case 2 % Step South
                    v1 = 0; v2 = -1;
                case 3 % Step East
                    v1 = 1; v2 = 0;
                case 4 % Step West
                    v1 = -1; v2 = 0;
                case 5 % No Movement
                    v1 = 0; v2 = 0;
            end
            old = var.i_nodes_coordinates(i,:);
            
            if(old(1) == 1 && v1 == -1)
                var.i_nodes_coordinates(i,1) = p.Nx;
            elseif(old(1) == p.Nx && v1 == 1)
                var.i_nodes_coordinates(i,1) = 1;
            else
                var.i_nodes_coordinates(i,1) = old(1) + v1;
            end
            
            if(old(2) == 1 && v2 == -1)
                var.i_nodes_coordinates(i,2) = p.Ny;
            elseif(old(2) == p.Ny && v2 == 1)
                var.i_nodes_coordinates(i,2) = 1;
            else
                var.i_nodes_coordinates(i,2) = old(2) + v2 ;
                
            end
            new = var.i_nodes_coordinates(i,:);
            var.i_nodes_array(old(1),old(2)) = var.i_nodes_array(old(1),old(2)) - 1;
            var.i_nodes_array(new(1),new(2)) = var.i_nodes_array(new(1),new(2)) + 1;
            
            
            
        end
        var.nodes_coordinates = [p.x(var.i_nodes_coordinates(:,1))', p.y(var.i_nodes_coordinates(:,2))'];
        
    case "Randoms"
        if(var.current_count == var.change_count)
            var.current_count = 0;
            N = var.num_nodes;
            Asize = [p.Nx,p.Ny];
            [Ir,Ic] = ind2sub(Asize,randperm(prod(Asize),N));
            var.i_nodes_coordinates = [Ir',Ic'];
            var.nodes_coordinates = [p.x(var.i_nodes_coordinates(:,1))', p.y(var.i_nodes_coordinates(:,2))'];
            
        else
            var.current_count = var.current_count + 1;
        end
    case {"Lagrangian" "Local"}
        var = update_Lagrange(var,p,psi_hat);
    case "Airplane"
        var.nodes_coordinates(:,1) = var.nodes_coordinates(:,1) + var.nodes_velocities(:,1);
        var.nodes_coordinates(:,2) = var.nodes_coordinates(:,2) + var.nodes_velocities(:,2);
        
        var.nodes_coordinates = mod((var.nodes_coordinates + pi),2*pi) - pi;
        
end
end

function [energy,rhs_hat]= compute_rhs_hat(t,psi_hat,p) %#ok<*INUSL>
% Compute the fft of the right-hand side of the equation.

% Dealias by setting the "middle" Fourier coefficients to zero.
% This is where Matlab stores the highest frequencies.
psi_hat(p.dealias_modes) = 0;
energy = 0;

% % Compute the nonlinear terms in physical space using dealiased versions.
% psi_x      = ifftn(              psi_x_hat,'symmetric'); % psi_x = -v
% psi_y      = ifftn(              psi_y_hat,'symmetric'); % psi_y =  u
% psi_x_lap  = ifftn(p.k_lap.*psi_x_hat,'symmetric');
% psi_y_lap  = ifftn(p.k_lap.*psi_y_hat,'symmetric');
% % rhs_hat = p.helm_inv.*p.k_lap_inv.*fft2(-psi_y.*psi_x_lap + psi_x.*psi_y_lap + p.f_hat);
% rhs_hat = p.k_lap_inv.*(fft2(-psi_y.*psi_x_lap + psi_x.*psi_y_lap) + p.f_hat);

% Basdevant formula (4 fft/ifft operations, 6 products). See:
%   P. Emami and J. C. Bowman,"On the Global Attractor of 2D Incompressible Turbulence with Random Forcing",
%   Journal of Differential Equations. 264, 4036-4066 (2018).
u1       =  ifftn(bsxfun(@times,p.iky,psi_hat),'symmetric'); % psi_y =  u
u2       = -ifftn(bsxfun(@times,p.ikx,psi_hat),'symmetric'); % psi_x = -v
rhs_hat = p.k_lap_inv.*(p.k1sq_m_k2sq.*fftn(u1.*u2) + p.k1k2.*fftn(u2.^2 - u1.^2) + p.f_hat);


end % =========== End function rhs_hat ============

function [u1,u2] = psi_converter(psi_hat, p)
% Converts stream function psi to velocity components u1 and u2
u1 =  ifftn(bsxfun(@times,p.iky,psi_hat),'symmetric'); % psi_y =  u
u2 = -ifftn(bsxfun(@times,p.ikx,psi_hat),'symmetric'); % psi_x = -v
end

function [psi] = velocity_converter(u1,u2,p)
% Converts velocity in physical space to stream function psi
%(given in Fourier space)
% curl(u1,u2) = -Lap psi
u1_hat = fft2(u1);
u2_hat = fft2(u2);
u2_x_hat = bsxfun(@times,p.ikx,u2_hat);
u1_y_hat = bsxfun(@times,p.iky,u1_hat);
psi =  -p.k_lap_inv.*(u2_x_hat - u1_y_hat);
end

function[x,y] = interpolated_points(var,p)
%Generates the points
switch(var.observer_type)
    case {"Lagrangian" "Airplane" "Local"}
        %         x_shift = ((p.Nx)/2+1);
        %         y_shift = ((p.Ny)/2+1); %????
        coor_x = (var.nodes_coordinates(:,1) + pi)./p.dx + 1;
        coor_y = (var.nodes_coordinates(:,2) + pi)./p.dy + 1;
        x1 = floor(coor_x);
        x2 = ceil(coor_x);
        y1 = floor(coor_y);
        y2 = ceil(coor_y);
        %         var.nodes_coordinates'
        %         [x1, x2, y1, y2]'
        x1 = var.indexing_array(x1);
        x2 = var.indexing_array(x2);
        y1 = var.indexing_array(y1);
        y2 = var.indexing_array(y2);
        x = [x1';x2'];
        y = [y1';y2'];
end

end
function[var] = update_Lagrange(var,p,psi_hat)
%Updates the Lagrangian trajectory
%l_n = dt l_{n-1} + dt(23/12 u_{n-1} - 4/3 u_{n-2} + 5/12 u_{n-3})
[u1,u2] = psi_converter(psi_hat,p);
u1_n = interp_from_grid(var,p,u1,var.nodes_coordinates);
u2_n = interp_from_grid(var,p,u2,var.nodes_coordinates);

var.nodes_coordinates(:,1) = var.nodes_coordinates(:,1) + var.velocity_x*p.dt*(23/12*u1_n + -4/3*var.u_old1_x + 5/12*var.u_old2_x);
var.nodes_coordinates(:,2) = var.nodes_coordinates(:,2) + var.velocity_y*p.dt*(23/12*u2_n + -4/3*var.u_old1_y + 5/12*var.u_old2_y);

var.nodes_coordinates = mod((var.nodes_coordinates + pi),2*pi) - pi;

var.u_old2_x =  var.u_old1_x;
var.u_old1_x = u1_n;
var.u_old2_y =  var.u_old1_y;
var.u_old1_y = u2_n;


end
function[z] = interp_from_grid(var,p,value,points)
%Interpolates off-grid function values from grid values
% Requires value (value to be interpolated) and (points) which are the
% query points
[x,y] = interpolated_points(var,p); % Grabs nearest points to the var
z_data = value(sub2ind(size(value),x,y)); %Generates z data for interpolant
x_pts = p.x(x)';
y_pts = p.y(y)';
F1 = scatteredInterpolant(x_pts,y_pts,z_data);
z = F1(points);
end
function[var] = reset(var)
var.resetTime = var.resetTime + var.resetInterval
switch var.resetConfig
    case 'random'
        var.nodes_coordinates = rand(var.num_nodes,2)*(2*pi) - pi;
        var.u_old1_x = zeros(length(var.nodes_coordinates),1);
        var.u_old2_x = zeros(length(var.nodes_coordinates),1);
        var.u_old1_y = zeros(length(var.nodes_coordinates),1);
        var.u_old2_y = zeros(length(var.nodes_coordinates),1);
end
end


function [u_obs] = compute_observations(u1, u2, var, p)
u_obs1 = zeros(p.Nx, p.Ny);
u_obs2 = zeros(p.Nx, p.Ny);
switch var.observer_type
    case {"Uniform", "Random Walkers", "Randoms"}
        x_pts_per = [];
        y_pts_per = [];
        u_data1 = [];
        u_data2 = [];
        r = var.i_nodes_coordinates(:,1);
        c = var.i_nodes_coordinates(:,2);
        x_pts = p.x(r)';
        y_pts = p.y(c)';
        
        obs_data_1 = u1(sub2ind(size(u1),r,c));
        obs_data_2 = u2(sub2ind(size(u2),r,c));
        %         u_obs1 = zeros(p.Nx, p.Ny);
        %         u_obs2 = zeros(p.Nx, p.Ny);
        
        for(i = -1:1)
            for(j=-1:1)
                x_pts_per = [x_pts_per; x_pts + p.Lx*i];
                y_pts_per = [y_pts_per; y_pts + p.Ly*j];
                
                u_data1 = [u_data1; obs_data_1];
                u_data2 = [u_data2; obs_data_2];
            end
        end
        
        
        F1 = scatteredInterpolant(x_pts_per,y_pts_per, u_data1);
        F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
        %         F1.Method = 'nearest';
        %         F2.Method = 'nearest';
        u_obs1 = F1(p.X,p.Y)';
        u_obs2 = F2(p.X,p.Y)';
        
        
        
        
    case "Thin Sweeps"
        
        u_obs1(:,var.i_nodes_y) = u1(:,var.i_nodes_y);
        u_obs2(:,var.i_nodes_y) = u2(:,var.i_nodes_y);
        
    case "Thin Sweeps (x2)"
        u_obs1 = zeros(size(u1));
        u_obs1(:,var.i_nodes_y) = u1(:,var.i_nodes_y);
        u_obs1(var.i_nodes_x,:) = u1(var.i_nodes_x,:);
        
        u_obs2 = zeros(size(u1));
        u_obs2(:,var.i_nodes_y) = u2(:,var.i_nodes_y);
        u_obs2(var.i_nodes_x,:) = u2(var.i_nodes_x,:);
    case {"Lagrangian" "Airplane"}
        x_pts_per = [];
        y_pts_per = [];
        u_data1 = [];
        u_data2 = [];
        obs_data_1 = interp_from_grid(var,p,u1,var.nodes_coordinates);
        obs_data_2 = interp_from_grid(var,p,u2,var.nodes_coordinates);
        x_pts = var.nodes_coordinates(:,1);
        y_pts = var.nodes_coordinates(:,2);
        
        for(i = -1:1)
            for(j=-1:1)
                x_pts_per = [x_pts_per; x_pts + p.Lx*i];
                y_pts_per = [y_pts_per; y_pts + p.Ly*j];
                
                u_data1 = [u_data1; obs_data_1];
                u_data2 = [u_data2; obs_data_2];
            end
        end
        
        
        F1 = scatteredInterpolant(x_pts_per,y_pts_per, u_data1);
        F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
        F1.Method = 'nearest';
        F2.Method = 'nearest';
        
        u_obs1 = F1(p.X,p.Y)';
        u_obs2 = F2(p.X,p.Y)';
        
    case "Thick Sweeps"
        x_pts_per = [];
        y_pts_per = [];
        u_data1 = [];
        u_data2 = [];
        r = var.i_nodes_coordinates(:,1);
        c = var.i_nodes_coordinates(:,2);
        x_pts = p.x(r)';
        y_pts = p.y(c)';
        
        x_zeros = p.x(var.zeros_coordinates(:,1))';
        y_zeros = p.y(var.zeros_coordinates(:,2))';
        zero_data_1 = zeros(length(x_zeros),1);
        zero_data_2 = zeros(length(x_zeros),1);
        
        obs_data_1 = u1(sub2ind(size(u1),r,c));
        obs_data_2 = u2(sub2ind(size(u2),r,c));
        
        x_pts = [x_pts; x_zeros];
        y_pts = [y_pts; y_zeros];
        obs_data_1 = [obs_data_1; zero_data_1];
        obs_data_2 = [obs_data_2; zero_data_2];
        
        for(i = -1:1)
            for(j=-1:1)
                x_pts_per = [x_pts_per; x_pts + p.Lx*i];
                y_pts_per = [y_pts_per; y_pts + p.Ly*j];
                
                u_data1 = [u_data1; obs_data_1];
                u_data2 = [u_data2; obs_data_2];
            end
        end
        
        
        F1 = scatteredInterpolant(x_pts_per,y_pts_per, u_data1);
        F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
        u_obs1 = F1(p.X,p.Y)';
        u_obs2 = F2(p.X,p.Y)';
        
    case "Local"
        x_pts_per = [];
        y_pts_per = [];
        u_data1 = [];
        u_data2 = [];
        obs_data_1 = interp_from_grid(var,p,u1,var.nodes_coordinates);
        obs_data_2 = interp_from_grid(var,p,u2,var.nodes_coordinates);
        x_pts = var.nodes_coordinates(:,1);
        y_pts = var.nodes_coordinates(:,2);
        
        for(i = -1:1)
            for(j=-1:1)
                x_pts_per = [x_pts_per; x_pts + p.Lx*i];
                y_pts_per = [y_pts_per; y_pts + p.Ly*j];
                
                u_data1 = [u_data1; obs_data_1];
                u_data2 = [u_data2; obs_data_2];
            end
        end
        
        F1 = scatteredInterpolant(x_pts_per,y_pts_per, u_data1);
        F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
        F1.Method = 'nearest';
        F2.Method = 'nearest';
        
        u_obs1(sub2ind(size(u_obs1),var.nearby_indices(:,1),var.nearby_indices(:,2))) = F1(var.nearby_points(:,1),var.nearby_points(:,2));
        u_obs2(sub2ind(size(u_obs2),var.nearby_indices(:,1),var.nearby_indices(:,2))) = F2(var.nearby_points(:,1),var.nearby_points(:,2));
        
end
u_obs = velocity_converter(u_obs1,u_obs2,p);

end

function [v_smoothing] = compute_Ih_v(var, p)
switch var.interp_type
    case "Standard"
        [v1,v2] = psi_converter(var.v_hat, p);
        v_smoothing = compute_observations(v1,v2, var, p);
end
end


function [vars] = initialize_default_var(load_vars, p)
% Checks var and initializes missing parameters to default values
%
vars = [];
vars = repelem(struct,1,1);
offset = 0;

for(i = 1: length(load_vars))
    var = load_vars(i);
    if(~isfield(var, 'observer_type'))
        var.observer_type = "Uniform";
        var.mu = 10;
        
        
        var.v = zeros([p.Nx,p.Ny]);
        var.v_hat = fft2(var.v);
        var.v2 = var.v_hat;
        var.v3 = var.v_hat;
        
        var.v1_mu = zeros([p.Nx,p.Ny]);
        
        
        var.int_nodes_x = 75; %75
        var.int_nodes_y = 75; %75
        var.i_nodes_x = floor(linspace(1,p.Nx,var.int_nodes_x));
        var.i_nodes_y = floor(linspace(1,p.Ny,var.int_nodes_y));
        
        [A,B] = meshgrid(var.i_nodes_x,var.i_nodes_y);
        c=cat(2,A',B');
        var.nodes_index =reshape(c,[],2);
        var.nodes_location = var.nodes_index;
        var.nodes_location(:,1) = p.x(var.nodes_index(:,1));
        var.nodes_location(:,2) = p.y(var.nodes_index(:,2));
        
        var.i_nodes_array = var.v;
        for j = 1: var.int_nodes_x
            for k = 1: var.int_nodes_y
                var.i_nodes_array(var.i_nodes_x(j),var.i_nodes_y(k)) = 1;
            end
        end
        var.i_nodes_array = sparse(var.i_nodes_array);
        [a,b,~] = find(var.i_nodes_array);
        var.i_nodes_coordinates = [a,b];
        var.nodes_coordinates = [p.x(var.i_nodes_coordinates(:,1))', p.y(var.i_nodes_coordinates(:,2))'];
        
        
        
        var.error = ones(1,p.nT_steps-offset)*nan;
        var.cpu_time = zeros(1,p.nT_steps);
        
        var.ens_u = ones(1,p.nT_steps-offset)*nan;
        var.ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
        var.ens_v = ones(1,p.nT_steps-offset)*nan;
        var.ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
        var.ens_umv = ones(1,p.nT_steps-offset)*nan;
        var.ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;
        
        
    end
    
    if(~isfield(var, 'interp_type'))
        var.interp_type = "Standard";
    end
    
    if(~isfield(var, 'time_type'))
        var.time_type = "None";
        var.obs_dt = 0;
        var.obs_time = 0;
        
    end
    
    if(~isfield(var, 'reset_type'))
        var.reset_type = "None";
        var.resetTimeInterval = NaN;
        var.resetTime = NaN;
    end
    
    if(~isfield(var, 'trunc_num'))
        var.trunc_num = inf;
        var.trunc_array = ones(p.Nx,p.Ny);
        
    end
    if(~isfield(var,'indexing_array'))
        var.indexing_array = [1:p.Nx,1:p.Nx];
        
    end
    if(~isfield(var,'v1_mu'))
        var.v1_mu = zeros(p.Nx,p.Ny);
        
    end
    %Load all fields into new struct
    %Have to do it this way to avoid dissimilar structures error.
    for( fn = fieldnames(var)')
        vars(i).(fn{1}) = var.(fn{1});
    end
    
end

end
