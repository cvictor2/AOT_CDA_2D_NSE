function vars = DataAssimilationVariables_NSE_CLEANED(p)
% Author: Collin Victor. Last modified on 2020-09-25.
% Adapted for use by NSE 2D DA programs
% Initializes information for various types of observer regimes
% Including: Movement types of observers (i.e. mobile and static), type of
% interpolation, time type of observers (i.e. Chi vs Delta), and observer reset
% functionality.



vars = [];
vars = repelem(struct,1,1);
i = 0;
offset = 0;

dx = p.dx;
dy = p.dy;
dt = p.dt;
Nx = p.Nx;
Ny = p.Ny;





% Uncomment the type of observer you want to use.
%% Standard uniform grid setup
      i = i+1;
    vars(i).observer_type = "Uniform";
    vars(i).interp_type = "Standard";
    vars(i).trunc_num = 100;
    vars(i).trunc_array = ones(Nx,Ny);
for sj = 1:Ny
    for si = 1:Nx
        if(p.kx(si)^2 + p.ky(sj)^2 >= vars(i).trunc_num^2)
            vars(i).trunc_array(si,sj) = 0;
        end
    end
end
    
    
%     vars(i).mu = 300;
   vars(i).mu = 15;
%     vars(i).mu = 1/3;

    
    
    
    vars(i).v = zeros([p.Nx,p.Ny]);
    vars(i).v_hat = fft2(vars(i).v);
    vars(i).v2 = vars(i).v_hat;
    vars(i).v3 = vars(i).v_hat;

    vars(i).v1_mu = zeros([p.Nx,p.Ny]);

    
    vars(i).int_nodes_x = 75;
    vars(i).int_nodes_y = 75;

    vars(i).i_nodes_x = floor(linspace(1,p.Nx,vars(i).int_nodes_x + 1));
    vars(i).i_nodes_y = floor(linspace(1,p.Ny,vars(i).int_nodes_y + 1));
    vars(i).i_nodes_x(end) = [];
    vars(i).i_nodes_y(end) = [];
    [A,B] = meshgrid(vars(i).i_nodes_x,vars(i).i_nodes_y);
    c=cat(2,A',B');
    vars(i).nodes_index =reshape(c,[],2);
    vars(i).nodes_location = vars(i).nodes_index;
    vars(i).nodes_location(:,1) = p.x(vars(i).nodes_index(:,1));
    vars(i).nodes_location(:,2) = p.y(vars(i).nodes_index(:,2));

    vars(i).i_nodes_array = vars(i).v;
    for j = 1: vars(i).int_nodes_x
        for k = 1: vars(i).int_nodes_y
            vars(i).i_nodes_array(vars(i).i_nodes_x(j),vars(i).i_nodes_y(k)) = 1;
        end
    end
    vars(i).i_nodes_array = sparse(vars(i).i_nodes_array);
    [a,b,~] = find(vars(i).i_nodes_array);
    vars(i).i_nodes_coordinates = [a,b];
    vars(i).nodes_coordinates = [p.x(vars(i).i_nodes_coordinates(:,1))', p.y(vars(i).i_nodes_coordinates(:,2))'];



    vars(i).error = ones(1,p.nT_steps-offset)*nan;
    vars(i).cpu_time = zeros(1,p.nT_steps);
    
    vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
    vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
    vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
    vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
    vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
    vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;

    
    vars(i).resetTimeInterval = NaN;
    vars(i).resetTime = vars(i).resetTimeInterval;
    vars(i).resetConfig = 'none';    

%% Sweep with only y
%         i = i+1;
%     vars(i).observer_type = "Thin Sweeps";
%     vars(i).interp_type = "Standard";
%     vars(i).v = zeros([p.Nx,p.Ny]);
%     vars(i).v_hat = fft2(vars(i).v);
%     vars(i).v2 = fft2(vars(i).v);
%     vars(i).v3 = vars(i).v2;
% 
% %     vars(i).i_nodes_y = [1:50];
%     vars(i).i_nodes_y = [1:10];
% %     vars(i).i_nodes_y = [1:19];
% %     vars(i).i_nodes_y = [1,150, 300];
%     vars(i).i_nodes_x = [1:1];
% 
% 
%     velocity_x  = (length(vars(i).i_nodes_x))*dx/dt;
%     velocity_y = 10*(length(vars(i).i_nodes_y))*dy/dt;
% %     velocity_y = (3)*dy/dt;
% 
%     vars(i).velocity_y = velocity_y;
%     vars(i).velocity_x = velocity_x;
%     vars(i).indexing_array = [1:Nx,1:1+ceil(max(velocity_x*dt/dx, velocity_y*dt/dy))];
%     
% %     vars(i).i_nodes_coordinates = ones(2*p.N+2*length(vars(i).i_nodes_y),2);
% %     vars(i).i_nodes_coordinates(1:p.N,1) = 1:p.N;
% %     vars(i).i_nodes_coordinates(1:p.N,2) = 1;
% %     vars(i).i_nodes_coordinates(p.N+1,1) = 1:p.N;
% %     vars(i).i_nodes_coordinates(2*p.N,2) = vars(i).i_nodes_y(end);
% %     vars(i).i_nodes_coordinates
% %     vars(i).nodes_coordinates = [p.x(vars(i).i_nodes_coordinates(:,1))', p.y(vars(i).i_nodes_coordinates(:,2))'];
% 
% 
%     vars(i).error = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;
% 
%     vars(i).mu = 50;
%     vars(i).cpu_time = zeros(1,p.nT_steps);

%
%% Type 3 Sweeping with X and Y
% % %     i = i+1;
% % %     vars(i).observer_type = "Thin Sweeps (x2)";
% % %     vars(i).v = zeros([p.Nx,p.Ny]);
% % % %     vars(i).int_nodes_x = 16;
% % % %     vars(i).int_nodes_y = 16;
% % % %     vars(i).i_nodes_x = floor(linspace(1,p.Nx-1,vars(i).int_nodes_x));
% % % %     vars(i).i_nodes_y = [1:50];
% % %     vars(i).i_nodes_y = unique([1:5:Ny]);
% % %     vars(i).i_nodes_x = [1:2];
% % % %     vars(i).i_nodes_x = unique( [1:10:Nx, Nx]);
% % %
% % %
% % %     velocity_x  = (length(vars(i).i_nodes_x))*dx/dt;
% % % %     velocity_y = (length(vars(i).i_nodes_y))*dy/dt;
% % %     velocity_y = 0;
% % %     vars(i).velocity_y = velocity_y;
% % %     vars(i).velocity_x = velocity_x;
% % %     vars(i).indexing_array = [1:Nx,1:1+ceil(max(velocity_x*dt/dx, velocity_y*dt/dy))];
% % %
% % %     vars(i).error = ones(1,p.nT_steps-offset)*nan;
% % %     vars(i).mu = 1000;
% %
%% Type 4: Random Walkers - Creeps
%         i = i+1;
% 
% 
%     vars(i).observer_type = "Random Walkers";
%     vars(i).v = zeros([p.Nx,p.Ny]);
%     vars(i).v_hat = fft2(vars(i).v);
%     vars(i).v2 = fft2(vars(i).v);
%     vars(i).v3 = vars(i).v2;
% %
%     vars(i).int_nodes_x = 5;
%     vars(i).int_nodes_y = 5;
% %     vars(i).int_nodes_x = 75;
% %     vars(i).int_nodes_y = 75;
% %     vars(i).int_nodes_x = 50;
% %     vars(i).int_nodes_y = 50;
% %     vars(i).int_nodes_x = 45;
% %     vars(i).int_nodes_y = 45;
% vars(i).num_nodes = vars(i).int_nodes_x*vars(i).int_nodes_y;
%     vars(i).i_nodes_x = floor(linspace(1,p.Nx,vars(i).int_nodes_x+1));
%     vars(i).i_nodes_y = floor(linspace(1,p.Ny,vars(i).int_nodes_y+1));
%     vars(i).i_nodes_x(end) = [];
%     vars(i).i_nodes_y(end) = [];
%     vars(i).i_nodes_array = vars(i).v;
%     for j = 1: vars(i).int_nodes_x
%         for k = 1: vars(i).int_nodes_y
%             vars(i).i_nodes_array(vars(i).i_nodes_x(j),vars(i).i_nodes_y(k)) = 1;
%         end
%     end
%     vars(i).i_nodes_array = sparse(vars(i).i_nodes_array);
%     [a,b,~] = find(vars(i).i_nodes_array);
%     vars(i).i_nodes_coordinates = [a,b];
%     vars(i).nodes_coordinates = [p.x(vars(i).i_nodes_coordinates(:,1))', p.y(vars(i).i_nodes_coordinates(:,2))'];
%     
%     velocity_x = 0;
%     velocity_y = 0;
%     vars(i).velocity_y = velocity_y;
%     vars(i).velocity_x = velocity_x;
%     vars(i).indexing_array = [1:Nx,1:Nx];
% 
%     vars(i).error = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;
% 
%     
%     
%     vars(i).mu = 10;
%     vars(i).cpu_time = zeros(1,p.nT_steps);
% 
%% Type 5: Randoms
%         i = i+1;
% 
% 
%     vars(i).observer_type = "Randoms";
%     vars(i).v = zeros([p.Nx,p.Ny]);
%     vars(i).v_hat = fft2(vars(i).v);
%     vars(i).v2 = fft2(vars(i).v);
%     vars(i).v3 = vars(i).v2;
% 
% %     vars(i).int_nodes_x = 35;
% %     vars(i).int_nodes_y = 35;
% 
% %     vars(i).int_nodes_x = 25;
% %     vars(i).int_nodes_y = 25;
%     vars(i).int_nodes_x = 75;
%     vars(i).int_nodes_y = 75;
% 
%     vars(i).num_nodes = vars(i).int_nodes_x*vars(i).int_nodes_y;
% 
%     N = vars(i).num_nodes;
%     Asize = [p.Nx,p.Ny];
%     [Ir,Ic] = ind2sub(Asize,randperm(prod(Asize),N));
%     vars(i).i_nodes_coordinates = [Ir',Ic'];
%     vars(i).nodes_coordinates = [p.x(vars(i).i_nodes_coordinates(:,1))', p.y(vars(i).i_nodes_coordinates(:,2))'];
%     
%     
%     
%     vars(i).error = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;
% 
%     
%     
%     vars(i).change_count = 1;
%     vars(i).current_count = 0;
%     vars(i).mu = 10;
%     vars(i).cpu_time = zeros(1,p.nT_steps);

% 
% 
%% Type 6: Lagrangian Particles
%         i = i+1;
% 
% 
%     vars(i).observer_type = "Lagrangian";
%     vars(i).interp_type = "Standard";
%     vars(i).v = zeros([p.Nx,p.Ny]);
%     vars(i).v_hat = fft2(vars(i).v);
%     vars(i).v2 = fft2(vars(i).v);
%     vars(i).v3 = vars(i).v2;
% 
% %     vars(i).int_nodes_x = 25;
% %     vars(i).int_nodes_y = 25;
%     
%     
% %     vars(i).int_nodes_x = 35;
% %     vars(i).int_nodes_y = 35;
% 
%     vars(i).int_nodes_x = 2;
%     vars(i).int_nodes_y = 5;
% %     vars(i).int_nodes_x = 100;
% %     vars(i).int_nodes_y = 100;
% 
%     vars(i).num_nodes = vars(i).int_nodes_x*vars(i).int_nodes_y;
% % %     vars(i).num_nodes = 2000;
% % %     vars(i).num_nodes = vars(i).int_nodes_x*vars(i).int_nodes_y;
% %     vars(i).i_nodes_x = floor(linspace(1,p.Nx,vars(i).int_nodes_x));
% %     vars(i).i_nodes_y = floor(linspace(1,p.Ny,vars(i).int_nodes_y));
% %     vars(i).i_nodes_array = vars(i).v;
% %     for j = 1: vars(i).int_nodes_x
% %         for k = 1: vars(i).int_nodes_y
% %             vars(i).i_nodes_array(vars(i).i_nodes_x(j),vars(i).i_nodes_y(k)) = 1;
% %         end
% %     end
% %     vars(i).i_nodes_array = sparse(vars(i).i_nodes_array);
% %     [a,b,~] = find(vars(i).i_nodes_array);
% %     vars(i).i_nodes_coordinates = [a,b];
% % %     vars(i).nodes_coordinates = [p.x(vars(i).i_nodes_coordinates(:,1))', p.y(vars(i).i_nodes_coordinates(:,2))'];
%     
%     vars(i).nodes_coordinates = rand(vars(i).num_nodes,2)*(2*pi) - pi;
% 
%     
%     
%     
% %     vars(i).num_nodes = 5625;
% %     vars(i).nodes_coordinates = rand(vars(i).num_nodes,2)*(2*pi) - pi;
%     
% 
%     vars(i).velocity_x = 1;
%     vars(i).velocity_y = 1;
%     vars(i).indexing_array = [1:Nx,1:Nx];
% 
%     vars(i).u_old1_x = zeros(length(vars(i).nodes_coordinates),1);
%     vars(i).u_old2_x = zeros(length(vars(i).nodes_coordinates),1);
%     vars(i).u_old1_y = zeros(length(vars(i).nodes_coordinates),1);
%     vars(i).u_old2_y = zeros(length(vars(i).nodes_coordinates),1);
% 
%     
%     vars(i).error = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;
%     
%     vars(i).mu = 10;
%     vars(i).cpu_time = zeros(1,p.nT_steps);
% 
%     vars(i).resetTimeInterval = inf;
%     vars(i).resetTime = vars(i).resetTimeInterval;
%     vars(i).resetConfig = 'random';

%% Type 7: Sweep with Uniform Box
% i = i+1;
% vars(i).observer_type = "Thick Sweeps";
% vars(i).v = zeros([p.Nx,p.Ny]);
% vars(i).v_hat = fft2(vars(i).v);
% vars(i).v2 = vars(i).v_hat;
% vars(i).v3 = vars(i).v2;
% 
% 
% vars(i).int_nodes_y = 5;
% vars(i).int_nodes_x = 20;
% % vars(i).int_nodes_y = 20;
% % vars(i).int_nodes_x = 60;
% % vars(i).int_nodes_y = 38;
% % vars(i).int_nodes_x = 150;
% %     vars(i).i_nodes_y = [1:50];
% % vars(i).i_nodes_y = [1,3,5];
% %     vars(i).i_nodes_y = [1:20];
% %     vars(i).i_nodes_y = [1,150, 300];
% vars(i).i_nodes_y = floor(linspace(1,p.Ny/4,vars(i).int_nodes_y));
% %     vars(i).i_nodes_x = [1:1];
% vars(i).i_nodes_x = floor(linspace(1,p.Nx,vars(i).int_nodes_x + 1));
% vars(i).i_nodes_x(end) = [];
% vars(i).num_nodes = vars(i).int_nodes_y*vars(i).int_nodes_x;
% 
% vars(i).zeros_y = [p.Ny, vars(i).i_nodes_y(end)+1];
% vars(i).zeros_x = 1:p.Nx;
% [A,B] = meshgrid(vars(i).zeros_x,vars(i).zeros_y);
% c=cat(2,A',B');
% vars(i).zeros_coordinates =reshape(c,[],2);
% 
% 
% 
% vars(i).i_nodes_array = vars(i).v;
% for j = 1: vars(i).int_nodes_x
%     for k = 1: vars(i).int_nodes_y
%         vars(i).i_nodes_array(vars(i).i_nodes_x(j),vars(i).i_nodes_y(k)) = 1;
%     end
% end
% vars(i).i_nodes_array = sparse(vars(i).i_nodes_array);
% [a,b,~] = find(vars(i).i_nodes_array);
% vars(i).i_nodes_coordinates = [a,b];
% vars(i).nodes_coordinates = [p.x(vars(i).i_nodes_coordinates(:,1))', p.y(vars(i).i_nodes_coordinates(:,2))'];
% 
% 
% velocity_x  = 0;
% % velocity_y = 1*(length(vars(i).i_nodes_y))*dy/dt;
% % velocity_y = 1*floor((p.Nx)*p.dt)*dy/dt;
% % velocity_y = 1*dy/dt;
% velocity_y = 1*dy/dt;
% 
% % velocity_y = (p.Nx/4)*dy/dt;
% 
% 
% vars(i).velocity_y = velocity_y;
% vars(i).velocity_x = velocity_x;
% vars(i).indexing_array = [1:Nx,1:1+ceil(max(velocity_x*dt/dx, velocity_y*dt/dy))];
% 
%     vars(i).error = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).mu = 30;
% vars(i).cpu_time = zeros(1,p.nT_steps);



%% Type 8: Airplane Particles
%         i = i+1;
% 
% 
%     vars(i).observer_type = "Airplane";
%     vars(i).v = zeros([p.Nx,p.Ny]);
%     vars(i).v_hat = fft2(vars(i).v);
%     vars(i).v2 = fft2(vars(i).v);
%     vars(i).v3 = vars(i).v2;
% 
% %     vars(i).int_nodes_x = 25;
% %     vars(i).int_nodes_y = 25;
% 
% %     vars(i).int_nodes_x = 10;
% %     vars(i).int_nodes_y = 10;
% %     vars(i).int_nodes_x = 100;
% %     vars(i).int_nodes_y = 100;
% 
% %     vars(i).num_nodes = vars(i).int_nodes_x*vars(i).int_nodes_y;
% %     vars(i).num_nodes = 5625;
% %     vars(i).num_nodes = 1000;
%     vars(i).num_nodes = 10;
% 
%     vars(i).nodes_coordinates = rand(vars(i).num_nodes,2)*(2*pi) - pi;
%     vars(i).nodes_velocities = (rand(vars(i).num_nodes,2)*2)-1;
%     
%     vars(i).velocity_x = 1;
%     vars(i).velocity_y = 1;
%     
%     vars(i).nodes_velocities(:,1) =  vars(i).nodes_velocities(:,1)*vars(i).velocity_x*p.dx;
%     vars(i).nodes_velocities(:,2) =  vars(i).nodes_velocities(:,2)*vars(i).velocity_y*p.dy;
% 
%     
%     vars(i).indexing_array = [1:Nx,1:Nx];
% 
% 
%     
%     vars(i).error = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_u_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_v_Linf = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv = ones(1,p.nT_steps-offset)*nan;
%     vars(i).ens_umv_Linf = ones(1,p.nT_steps-offset)*nan;    vars(i).change_count = 1;
%     vars(i).current_count = 0;
%     vars(i).mu = 10;
%     vars(i).cpu_time = zeros(1,p.nT_steps);
%  

end