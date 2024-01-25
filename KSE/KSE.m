function KSE()
close all;
N = 2^8;
Lx = 32*pi;
dx = Lx/N;
x = 0:dx:Lx - dx;

lambda = 1;

T = 30;
dt = 0.01;
% dt = 1.2207e-4;
show = 10;

mu = 100;

k = [0:N/2-1 0 -N/2+1:-1]*(2*pi/Lx);
E = exp(dt*(lambda*k.^2 - k.^4));

dealias_mask = abs(k) <= floor((2/3)*N);
% [~,dealias_modes] = find(abs(k) > floor(N*2/3));


num_timesteps = ceil(T/dt);
t = dt:dt:dt*num_timesteps;

u_hat = zeros(N,1);

u_0 = cos(x/16).*(1+sin(x/16));

u_hat = fft(u_0);


spec = generate_spectrum_1D(u_hat);
modes = 1:N/2;


soln_history = zeros(N, ceil(num_timesteps/show)+1);
soln_history(:,1) = u_0;
plot_time = [0];

ref_fig = figure;
ref_soln_plot = scatter(x, ifft(u_hat,'symmetric'));
axis([0, Lx, -3,3]);

spec_fig = figure;
spec_plot = loglog(modes, spec);

error_fig = figure;

var = [];
observed_modes = 20;

trunc_array = zeros(N,1);

for i = 1:N
    if(abs(k(i)./(2*pi/Lx)) < observed_modes)
        trunc_array(i) = 1;
    end

end
trunc_array(1) = 0;
trunc_array(N/2+1) = 0;

trunc_index = find(trunc_array == 1);
trunc_index_comp = find(trunc_array == 0);

ensemble_count = 10000;
variance = 0.0000000001;



observations = u_hat(trunc_index);






ramp_up_timesteps = floor(50/dt);

error = NaN(1,num_timesteps);
error_aot = error;

for ti = 1:ramp_up_timesteps

    % u_dealiased = ifft(u_hat.*dealias_mask,'symmetric');
    % u_dealiased = ifft(u_hat,'symmetric');

    nonlin_term = (1i*k/2).*fft(real(ifft(u_hat.*dealias_mask)).^2);
    % nonlin_term = fft(u_dealiased.*real(ifft(1i*k.*u_hat, 'symmetric')));

    u_hat = E.*(u_hat - dt*nonlin_term);

end









enable_EnKF = true;



M = ensemble_count;
ensemble = struct('noise',[],'forecast',[],'analysis',[],'obs',[]);
ensembles = repmat(ensemble,1,M);
for j = 1:M
    ensembles(j) = initialize_ensemble(variance, N, trunc_index, trunc_index_comp,u_hat);

end

% figure;
% for j = 1:M

    % ens_phys = real(ifft(ensembles(j).forecast,'symmetric'));
    % plot(x, ens_phys);
    % hold on;

% end



tol = .1;
reset_times = [];

% AOT (nudging) solution
aot_hat = zeros(size(u_hat));
aot_hat(trunc_index) = u_hat(trunc_index);

n = 1;
for ti = 1:num_timesteps
    u_hat_old = u_hat;
    % u_dealiased = ifft(u_hat.*dealias_mask,'symmetric');
    % u_dealiased = ifft(u_hat,'symmetric');
      nonlin_term = (1i*k/2).*fft(real(ifft(u_hat.*dealias_mask)).^2);
    % nonlin_term = fft(u_dealiased.*real(ifft(1i*k.*u_hat, 'symmetric')));

    u_hat = E.*(u_hat - dt*nonlin_term);
  

    %observe previous timestep for nudging
    aot_obs = u_hat_old;
    %Zero out unobserved modes on observation data
    aot_obs(trunc_index_comp) = 0;

    % aot_obs(trunc_index_comp) = 0;

    %compute nudging feedback term I_h(u-v)
    Ihumv = aot_obs - aot_hat;
    %Zero out all unobserved modes
    Ihumv(trunc_index_comp) = 0;
    


    nonlin_aot = (1i*k/2).*fft(real(ifft(aot_hat.*dealias_mask)).^2);
    aot_hat = E.*(aot_hat - dt*nonlin_aot + dt*mu*(Ihumv));

    obs = u_hat.';


    if enable_EnKF
    %% Apply the EnKF
    total_forecast = zeros(N,1);
    total_noise = zeros(2*observed_modes - 2,1);
    total_analysis = zeros(N,1);
    for j = 1:M


            error_ens = norm(abs(ensembles(j).obs(trunc_index) - ensembles(j).analysis(trunc_index)),'fro')/N;
            if mod(t(ti),1) <= dt && error_ens >= tol
            % if  error_ens >= tol

                ensembles(j) = initialize_ensemble(variance, N, trunc_index, trunc_index_comp,u_hat);
                % reset = reset + 1;
                reset_times = unique([reset_times;t(ti)]);
            end

        %Draw ensemble observations
        fake_noise = generate_noise(variance,N);
        ensembles(j).noise = fake_noise(trunc_index);
        ensembles(j).obs = zeros(N,1);
        ensembles(j).obs(trunc_index) = obs(trunc_index) + fake_noise(trunc_index);

        total_forecast = total_forecast + ensembles(j).forecast;
        total_noise = total_noise + ensembles(j).noise;

    end
    mean_noise = total_noise./M;
    mean_forecast = total_forecast./M;

    X = zeros(length(trunc_index),M);
    % X = zeros(N*p.Ny, M);
    Y = zeros(length(trunc_index),M);
    %Construct matrices
    for j = 1:M
        X(:,j) = (ensembles(j).forecast(trunc_index) - mean_forecast(trunc_index)).';
        Y(:,j) = X(:,j) - ((ensembles(j).noise - mean_noise));
        % X(:,j) = vars(i).ensemble(j).forecast(:) - mean_forecast(:);
        % Y(:,j) = X(vars(i).trunc_index,j) - (vars(i).ensemble(j).noise - mean_noise);
    end
    X = X;
    Y = Y;
    X = X./(sqrt(M-1));
    Y = Y./(sqrt(M-1));
    obs_count = sum(trunc_array(:));

    A = Y*Y.';
    B = X*Y.';
    K = B/A;
    % condest(K)

    %Update the ensemble
    for j =1:M
        %x^a = x^f + K(y-H(x^f))
        %Here x^f is a full matrix, and K only affects observed
        %modes
        ensembles(j).analysis = ensembles(j).forecast;
        % corrector_flat = K*(vars(i).ensemble(j).obs - vars(i).ensemble(j).forecast(vars(i).trunc_index));
        % corrector = reshape(corrector_flat, N, p.Ny);
        % vars(i).ensemble(j).analysis = vars(i).ensemble(j).forecast + corrector;
        ensembles(j).analysis(trunc_index) = ensembles(j).analysis(trunc_index) + K*(ensembles(j).obs(trunc_index) - ensembles(j).forecast(trunc_index));

        total_analysis = total_analysis + ensembles(j).analysis;



        %Evolve forecast forward in time to next timestep
        % [~,vars(i).ensemble(j).forecast] = evolve(p, ev, vars(i).ensemble(j).analysis);
        v_hat = ensembles(j).analysis.';
        nonlin_term = (1i*k/2).*fft(real(ifft(v_hat.*dealias_mask)).^2);
        % nonlin_term = fft(u_dealiased.*real(ifft(1i*k.*u_hat, 'symmetric')));

        v_hat = E.*(v_hat - dt*nonlin_term);
        ensembles(j).forecast = v_hat.';
    end
    error(ti) = norm(abs(obs- mean_forecast),'fro')/N;
    end


    error_aot(ti) = norm(abs(u_hat - aot_hat),'fro')/N;


 





    if mod(ti,show)==0
        % u_phys = real(ifft(ensemble(1).forecast,'symmetric'));
        u_phys = real(ifft(u_hat,'symmetric'));

        soln_history(:,n+1) = u_phys;
        plot_time = [plot_time; ti*dt];
        n = n+1;

        spec = generate_spectrum_1D(u_hat);
        figure(spec_fig);
        spec_plot = loglog(modes, spec);
        hold on

        spec_aot = generate_spectrum_1D(aot_hat);
        loglog(modes, spec_aot);
        
        max(abs(spec - spec_aot))
        
        if enable_EnKF
        spec_EnKF = generate_spectrum_1D(mean_forecast);
        
        loglog(modes, spec_EnKF);
        end
        hold off;


        title(sprintf('Energy spectrum at t = %1.2f',t(ti)));


        figure(ref_fig);
        hold off;
        ref_soln_plot = plot(x, ifft(u_hat,'symmetric'));
        hold on;

        aot_soln_plot = plot(x, ifft(aot_hat,'symmetric'));
        % plot(x, ifft(abs(aot_hat - u_hat), 'symmetric'));
        
        hold off;

        title(sprintf('Reference solution at t = %1.2f',t(ti)));
        axis([0, Lx, -3,3]);
        
        if enable_EnKF
        hold on;
        plot(x, ifft(mean_forecast,'symmetric'));
        hold off;
        end
        % drawnow;

        figure(error_fig);
        semilogy(dt:dt:dt*length(error), error);
        hold on;
        semilogy(dt:dt:dt*length(error_aot), error_aot);
        hold off;
        title("Error computed over time");
        drawnow;
    end

end

if(plot_time(end)~= t(end))
    u_phys = real(ifft(u_hat,'symmetric'));

    soln_history(:,end) = u_phys;
    plot_time = [plot_time; t(end)];


end


figure;
% surf(soln_history);
% surf(soln_history);
% shading interp;
% surf([0,t],x,soln_history), shading interp, lighting phong, axis tight
surf(plot_time,x,soln_history), shading interp, lighting phong, axis tight
view([-90 90]), colormap(autumn);
% colormap(autumn);
light('color',[1 1 0],'position',[-1,2,2])
material([0.30 0.60 0.60 40.00 1.00]);
end

function spectrum = generate_spectrum_1D(soln_hat)

[N] = length(soln_hat);
spectrum = zeros(1, N/2);


for j = 1:N/2
    % for i = 1:Nx/2
    modes = floor(abs(j));
    if modes <= N/2
        spectrum(modes) = spectrum(modes) + abs(soln_hat(j))^2;
    end
    % end
end
modes = 1:N/2;
spectrum = sqrt(spectrum)/N^2;

spectrum = max(spectrum, eps);
end

function noisy_obs = generate_noise(variance,N)
% % psi_hat_size = size(psi_hat);
% % phi_v1 = p.noise_level*randn(psi_hat_size);
% % phi_v2 = p.noise_level*randn(psi_hat_size);
%
% % phi_v1_hat = fft2(phi_v1);
% % phi_v2_hat = fft2(phi_v2);
%
% noise = p.noise_level * (randn([N,p.Ny]));% + 1i*randn(size(p.ikx)));
%
% noise = fft2(noise);
%
% % Generate Gaussian white noise for potential field in Fourier space
% % noise = p.noise_level * (randn([N,p.Ny]));% + 1i*randn(size(p.ikx)));
% noise(1,1) = 0;
% noise(p.trunc_array == 1) = 0;

% Define the size of your 1D domain
% N = p.Nx;



% Generate the noise in Fourier space
noise_fourier_space = variance*(randn(N/2,1)+1i*(randn(N/2, 1)));
% noise_fourier_space = settings.noise_level*(randn(N, N/2+1) + 1i * randn(N, N/2+1));

% Set the 0 wave mode to 0
noise_fourier_space(1) = 0 + 0i;

% Create the complex conjugate mirrored version of your noise
noise_fourier_space_full = zeros(N, 1);
% noise_fourier_space_full(:, 1:N/2+1) = noise_fourier_space;
% Flipud function flips the array up down, and fliplr function flips it left to right


noise_fourier_space_full(1:N/2) = noise_fourier_space;
pos_freq = 2:N/2;

noise_fourier_space_full(N - pos_freq + 2) =  conj(noise_fourier_space(pos_freq));
% conj function takes the complex conjugate
% noise_fourier_space_full(:, N/2+2:N) = rot90(conj(noise_fourier_space(:, 2:N/2)),2);

% %Rescale noise based on wave mode (pink noise)
% % Define the wave numbers kx and ky
% % [kx, ky] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
% % k = sqrt(kx.^2 + ky.^2);
% k = p.kx.^2 + p.ky.^2;
% % Define your scaling function and apply it to the noise
% scaling = (k ~= 0) ./ (k + (k==0));  % avoid division by zero at (0,0), this makes the scaling to 1/k
% noise_fourier_space_full = noise_fourier_space_full .* scaling;
%
% % Inverse Fourier transform the noise back to real space
% noise_real_space = ifft2(noise_fourier_space_full, 'symmetric');

noisy_obs = noise_fourier_space_full;


end

function ensemble = initialize_ensemble(variance, N, trunc_index, trunc_index_comp,u_hat)
    noise = generate_noise(variance,N);
    ensemble = struct();
    ensemble.noise = noise(trunc_index);
    ensemble.forecast = u_hat.';
    ensemble.forecast(trunc_index_comp) = 0;
    ensemble.forecast(trunc_index)= ensemble.forecast(trunc_index) + noise(trunc_index);
    ensemble.analysis = zeros(N,1);
    ensemble.obs = u_hat;
    ensemble.obs(trunc_index_comp) = 0;
end

