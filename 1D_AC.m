%function Chafee_Infante(N1,T,nui,alphai,graph)
%Allen-Cahn eq. u_t - nu*u_xx = alpha u - beta*u^3
%basic constants defined by problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solves Allen-Cahn equation using data assimilation on a uniform grid
%
%Meant to use loaded data from Data directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all;
close all;
clear all;
% param = struct; %All of the parameters for this problem: alpha, nu, min_nodes, L, N, mu

%% Load Variables
% load('Data\2018.01.23.01.24.03.575.mat')
% load('Data/Overhauled/2018.04.10.23.18.30_Initial_Ramp_Only_nu=1.0e-04_alpha=1_mu=100.mat')
graph = true;           %Display Graphs
% graph2 = false;           %Display Graphs
save_check = false;     %Save data after ramp up and at end
reset_vars = true;      %Override loaded variable configurations with new definitions
if( exist('param','var') == 0)   
    alpha = 1; % alpha is in [0.0001, 150]
    beta = 1;
    nu = 7.5e-6;
%     nu = 5e-5;
%     nu = 1e-4;
%     nu = 0.01;
    %     nu is in [7.5e-6,0.01]
    
    L = 1;
    N = 2^12;
%     mu = 451.5;
%     mu = 3560;
    mu = 500;
%     (2e3)-5.875;
    %mu in (420, 425]
    %min mu = 425, if we don't consider late structure collapse, 460 otherwise.
    min_nodes = max(ceil(.25*nu^-.5),5);
    
%     seed = randi(10000)
%     seed = 3372
%     seed = 5938
%     seed = 9183
%     seed = 1973 %seed for local error plot of stationary probe
%     seed = 7620 %seed for local error plot of stationary probe
    seed = 4368
% seed = 6898 % structure collapse at t = 48
%     seed = 5425
%          seed = 3552
         T_ramp = 50;
%     T = 500000;
     T =  25;
    %     param = struct('alpha',alpha,'nu',nu,'L',L,'mu',mu,'min_nodes',min_nodes,'seed',seed,'T',T);
    
    
    %% Initialize other constants
    
    dx = L/N;
    %x = linspace(0,L,N);
    x = 0 + dx*(0:(N-1));
    t_ramp = 0;
    
%     dt = .2;
dt = 0.001;
%0.001;
%     speed = 10/2;
    timesteps = ceil(T/dt);
    
    
    %% Initial Conditions
    % seed = randi(10000);
    % fprintf('seed: %d\n', seed);
    rng(seed);
    a = (10^(0))*randn(1,floor(N/5));
    k = 1:1:floor(N/5);
    u_0 = zeros(1,N);
%     u_0 = randn(1,N);
    for i = 1:floor(N/5)
        u_0 = u_0 + (a(i)*sin(2*pi*k(i)*x/L));
    end
% u_0 = 5*x;
    u_0 = u_0/(norm(u_0)*sqrt(dx))*10^-3;
    
    % dt = dx/(1/sqrt(alpha)); % CFL, Courant Fredrich Lewy number for advection

    
    
    u = u_0';
    ui = u(2:N);
    t=0;
    
    %% Initialize Param
    param = struct('alpha',alpha,'beta',beta,'nu',nu,'L',L, 'N',N,'mu',mu,'min_nodes',min_nodes,'seed',seed,'T',T, 't', t, 'dt', dt, 'x', x, 'dx', dx, 'u_0', u_0, 't_ramp',t_ramp,'t_adj',t-t_ramp,'u_ramp',u_0, 'timesteps',timesteps);
    
    %% Ramp up
    offset = 0;
    u_solo = struct('u',u,'type',0);
    theta = 2*alpha^2/beta^2;
    % Convex Splitting Method
%     A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (1+2*alpha)*speye(N-1);
%     A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (beta+2*alpha^2/beta^2)*speye(N-1);
%     B = speye(N-1) - dt*A;
%    B = I -dt\nu D_{xx} + \beta+2\alpha/\beta I
%   Scheme is BU^{n+1} = U^n
%     B = speye(N-1) + dt*(theta-1)*beta*speye(N-1) - dt*(nu/(dx^2))*gallery('tridiag',N-1,1,-2,1);
    B = speye(N-1) + dt.*((theta-1).*alpha).*speye(N-1) - dt*(nu/(dx^2))*gallery('tridiag',N-1,1,-2,1); %FIX ME, should not be 0.
    B = sparse(B);
%     C = inv(B);
% D =  (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1);

    H = ones(N-1,1);

    while(max(abs(u))<.8/sqrt(alpha/beta)&&t<T_ramp-dt)
%     while(true)
        if (graph && mod(offset,100) == 0)
            if(~Graphing_Overhaul(param,u_solo,[]))
                return
            end
            %     pause;
        end
        
        %solves pde
        offset = offset +1;
        t = t+dt;
        param.t = t;
%         param.t_adj = t - t_ramp;

        %     umv = u - v;
        %error_DA(k) = sqrt(dx)*norm(umv);
        ui = B\(ui.*(H+dt.*( -beta.*((ui.^2)) +(theta)*alpha.*H)) ); %Convex Splitting Method %%FIX ME based on above change to B.
        u(2:N) = ui;
        u_solo.u = u;
%         u_f = abs(fft(u)/N);
        %     umv_error(offset) = sqrt(dx)*norm(u-v);
    end
    t_ramp = t;
    param.t_ramp = t_ramp;
    param.t_adj = t-t_ramp;
    u_ramp = u;
    param.u_ramp = u_ramp;
    
    if(save_check)
        date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
        save(sprintf('Data/Overhauled/%s_Initial_Ramp_Only_nu=%.1d_alpha=%.0d_mu=%.0d.mat',date_string,nu,alpha,mu),'param');
    end
    vars = [];
else
    %% Restore from loaded variables
    alpha = param.alpha;
    beta = param.beta;
    nu = param.nu;
    L = param.L;
    N = param.N;
    mu = param.mu;
    min_nodes = param.min_nodes;
    
    seed = param.seed;
    T = param.T;
    t = param.t;
    
    x = param.x;
    dx = param.dx
    t_ramp = param.t_ramp;
    
    rng(seed);
    
    u_0 = param.u_0;
    u_ramp = param.u_ramp;
    u = u_ramp;
    u_solo = struct('u',u_ramp,'type',0);
    
    dt = param.dt;
    timesteps = param.timesteps;
    
    offset = t_ramp/dt;

    if(graph)
        if(~Graphing_Overhaul(param,u_solo,[]))
            return
        end
    end
    
    % Convex Splitting Method
%     A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (1+2*alpha)*speye(N-1);
%     A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (1+2*alpha^2/beta^2)*speye(N-1);

%     B = speye(N-1) - dt*A;
    
end
%% Set up variables to track
offset = 0; %comment this out is you only want to run to time T, not time t_ramp + T
if(( exist('vars','var') == 0)||reset_vars)
    %Type Descriptions:
    %   0 - Normal solution, no data assimilation
    %   1 - Uniform Grid Data Assimilation
    %   2 - Standard Car Data Assimilation
    %  -1 - Retro Car, old car method
    %   3 - Car with uniform grid along length (Hybrid Method)
    vars = [];
    vars = repelem(struct,1,1);
    i = 0;
% 
%         Standard uniform grid setup
    i = i+1;
    vars(i).type = 1;
    vars(i).v = zeros(size(u));
    vars(i).int_nodes = 100;
%     vars(i).int_nodes = 90;
    vars(i).i_nodes = floor(linspace(1,N-1,vars(i).int_nodes));
    vars(i).x_nodes = x(vars(i).i_nodes);
    vars(i).error = zeros(1,timesteps-offset+1);
    vars(i).gamma = 0.001;
%     
%             i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
% %     vars(i).int_nodes =20;
% %     vars(i).int_nodes =N/2;
%     vars(i).int_nodes =30;
% 
% 
% %     vars(i).i_nodes = 1:vars(i).int_nodes;
% %     vars(i).i_nodes = 1900:(1900-1)+vars(i).int_nodes;
%     vars(i).i_nodes = N/16:(N/16-1)+vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
% %     vars(i).velocity =1/8*dx/dt;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);

%     i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
% %     vars(i).int_nodes =20;
% %     vars(i).int_nodes =N/2;
%     vars(i).int_nodes =10;
% 
% 
% %     vars(i).i_nodes = 1:vars(i).int_nodes;
% %     vars(i).i_nodes = 1900:(1900-1)+vars(i).int_nodes;
%     vars(i).i_nodes = N/16:(N/16-1)+vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
% %     vars(i).velocity =1/8*dx/dt;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
%     vars(i).gamma = 0.03;



    
%      i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
% %     vars(i).int_nodes =20;
% %     vars(i).int_nodes =N/2;
%     vars(i).int_nodes =100;
% 
% 
% %     vars(i).i_nodes = 1:vars(i).int_nodes;
% %     vars(i).i_nodes = 1900:(1900-1)+vars(i).int_nodes;
%     vars(i).i_nodes = N/16:(N/16-1)+vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
% %     vars(i).velocity =1/8*dx/dt;
%     vars(i).velocity = 100*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
%     vars(i).gamma = 0.3;
% 

    
                i = i+1;
    vars(i).type = 4;
    vars(i).v = zeros(size(u));
%     vars(i).int_nodes =20;
%     vars(i).int_nodes =N/2;
    vars(i).int_nodes =100;


%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).i_nodes = 1900:(1900-1)+vars(i).int_nodes;
    vars(i).i_nodes = N/16:(N/16-1)+vars(i).int_nodes;
    vars(i).x_nodes = x(vars(i).i_nodes);
    vars(i).x_nodes_o = vars(i).x_nodes;
%     vars(i).velocity =1/8*dx/dt;
    vars(i).velocity = 100*dx/dt;
    vars(i).direction = 1;
    vars(i).error = zeros(1,timesteps-offset+1);

    vars(i).gamma = .01;
    
    
    
%     
%     
% %         Standard car setup
%     i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes =5;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
%         i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes =10;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
%     i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes =25;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
%     i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes =50;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
%         i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes =100;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).x_nodes_o = vars(i).x_nodes;
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
% 
% 


    
    %Templates (Don't adjust these, just make copies if you want to use an
    %           instance)
%     % Standard uniform grid setup
%     i = i+1;
%     vars(i).type = 1;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes = min_nodes;
%     vars(i).i_nodes = floor(linspace(1,N-1,vars(i).int_nodes));
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).error = zeros(1,timesteps-offset+1);
%     
%     % Standard car setup
%     i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes = 30;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).velocity = 10*dx/dt;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
%     %Retro Car Configuration
%     i = i+1;
%     vars(i).type = -1;
%     vars(i).v = zeros(size(u));
%     vars(i).length_car = 10;
%     vars(i).int_nodes = vars(i).length_car;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).velocity = vars(i).length_car;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
    
%     %Uniform Grid on Car Configuration
%     i = i+1;
%     vars(i).type = 3;
%     vars(i).v = zeros(size(u));
%     vars(i).length_car = 10; %Length of car (in gridpoints)
%                              %Set to N to cover entire grid
%     vars(i).int_nodes = 2;
%     vars(i).i_nodes = floor(linspace(1,vars(i).length_car, vars(i).int_nodes));
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).velocity = 10;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
%     clear i;
end

% P = parallel.pool.Constant(param);
% P.Value.
%% Solve System
for j = 1:timesteps-offset
    t = t+dt;
    param.t = t;
    param.t_adj = param.t_adj + dt;
%     U = parallel.pool.Constant(u);
%     U.Value;
    %update the v's
    
%     gamma = 0.1;
    indexingArray = [1:N-1,1:N-1];
    if(i~=0)
        for (k = 1:length(vars))
            indices = vars(k).i_nodes;
%             vars(k).error(j) = sqrt(dx)*norm(u-vars(k).v);
            vars(k).error(j) = sqrt(dx)*norm(u(indices)-vars(k).v(indices));
%             error = vars(k).error(j)
            switch vars(k).type
                case {2,-1, 3,4}
                    %update car nodes
                    vars(k) = Update_Overhaul(param, vars(k));
            end
            umv = u - vars(k).v;
            Ih_umv = LinTerp(param, vars(k), umv);
            if(vars(k).type == 4 || vars(k).type == 1)
                Ih_umv = Ih_umv/((sqrt(dx)*norm(Ih_umv))^vars(k).gamma);
%                 for(val = 1:length(vars(k).i_nodes))
% %                 for(val = 1:length(Ih_umv))
%                     if(Ih_umv(indexingArray(vars(k).i_nodes(val)))~=0)
%                         Ih_umv(indexingArray(vars(k).i_nodes(val))) = Ih_umv(indexingArray(vars(k).i_nodes(val)))/(abs(Ih_umv(indexingArray(vars(k).i_nodes(val))))^vars(k).gamma);
%                     end
%                 end
            end
%             Ih_umv(indexingArray(vars(k).i_nodes)) = Ih_umv(indexingArray(vars(k).i_nodes))./(abs(Ih_umv(indexingArray(vars(k).i_nodes))).^gamma);
            
%                     ui = C*(ui.*(H+dt.*( -beta.*((ui.^2)) + alpha.*H)) ); %Convex Splitting Method

            vars(k).v(2:N) = B\(vars(k).v(2:N)+dt*(-beta*vars(k).v(2:N).^3 + (theta)*alpha*(vars(k).v(2:N)) + (mu)*Ih_umv)); %Convex Splitting Method
%             B\(ui.*(ones(N-1,1)+dt.*(alpha.*((ui.^2)) + (theta).*beta.*ones(N-1,1) )));
        end
        
%    u(2:N) = B\((u(2:N).*(ones(N-1,1)+dt.*(-beta.*((u(2:N).^2)) + theta*alpha.*ones(N-1,1) ))));%Convex Splitting Method
               
  
    end
    %update u
%     ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(1-2*alpha^2/beta^2-(ui.^2)))));

    u(2:N) = B\((u(2:N).*(ones(N-1,1)+dt.*(-beta.*((u(2:N).^2)) + theta*alpha.*ones(N-1,1) ))));%Convex Splitting Method
    u_solo.u = u;
    
%     show = 100;
show = 100;
    if(j>2000)
        show = 250;
    end
    %     if(max(u)> 0.9/sqrt(alpha)&&graph && mod(k,5) == 0
    if(graph && mod(j,show) == 0)
        if(i~=0)
            if(~Graphing_Overhaul(param,u_solo,vars))
                return
            end
        else
            if(~Graphing_Overhaul(param,u_solo,[]))
                return
            end
        end
    end
end
%Get last error measurement
if(i~=0)
 for k = 1:length(vars)
      indices = vars(k).i_nodes;
            vars(k).error(j) = sqrt(dx)*norm(u-vars(k).v);
%       vars(k).error(end) = sqrt(dx)*norm(u(indices)-vars(k).v(indices));
        vars(k).error(end) = sqrt(dx)*norm(u-vars(k).v);
 end
end
 
 %% Plot Error
 figure
 legendInfo = cell(1,length(vars));
 if(i~=0)
 for k = 1:length(vars)
     if(k==2)
         semilogy(0:dt:T,vars(k).error,'--','LineWidth',1);
     else
         semilogy(0:dt:T,vars(k).error,'LineWidth',1);
     end
     hold on;
     switch vars(k).type
         case 0
             text = 'u';
         case 1
             text = 'Uniform Grid';
         case 2 
             text = 'Sweeping Probe';
         otherwise
             text = 'Misc.';
     end
     legendInfo{k} = [text];
 end
 
 %%%%%%%%%%%%%%%%%%% TEST
 % Find first time when error drops below 1e-8 (should do 1e-15)
%  for(i=1:length(vars))
%      disp(dt*find(vars(i).error<1e-8,1));
%  end
 end
 %%%%%%%%%%%%%%%%%%% TEST
legend(legendInfo);
xlabel('Time');
ylabel('L^2 Norm of u-v');
title('Error Plot for Different Methods');
% date = datestr(now,'yyyy.mm.dd.HH.MM.SS');
if(save_check)
    if(~exists(date_string))
        date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
    end
    save(sprintf('Data/Overhauled/%s_Initial_Ramp_Var_Config_nu=%.1d_alpha=%.0d_mu=%.0d.mat',date_string,nu,alpha,mu),'param','vars');
end

function [lin] = LinTerp(p,v,umv)
if(length(v.i_nodes)>1)
    switch v.type
        case 1
            lin = interp1(v.x_nodes,umv(v.i_nodes),p.x)';
            lin(end) = [];
        case {2,-1, 3, 4}
            
            [~,index] = max(abs(v.x_nodes(2:end)-v.x_nodes(1:end-1)));
            %         if(maximum*p.dx) > ceil((p.L/v.int_nodes)*p.dx)
            if(index~=1)
                group1 = struct('i_nodes',v.i_nodes(1:index),'x_nodes',v.x_nodes(1:index));
                lend1 = (group1.i_nodes(1)-2)*p.dx;
                rend1 = (group1.i_nodes(end)+2)*p.dx;
                group2 = struct('i_nodes',v.i_nodes(index+1:end),'x_nodes',v.x_nodes(index+1:end));
                lend2 = (group2.i_nodes(1)-2)*p.dx;
                rend2 = (group2.i_nodes(end)+2)*p.dx;
                
                if(group1.i_nodes(1)==1)
                    lin1 = interp1([group1.x_nodes,rend1],[umv(group1.i_nodes)',0],p.x,'linear',0)';
                else
                    lin1 = interp1([0,group1.x_nodes,rend1],[0,umv(group1.i_nodes)',0],p.x,'linear',0)';
                end
                if(group2.i_nodes(end)==p.N)
                    lin2 = interp1([lend2,group2.x_nodes],[0,umv(group2.i_nodes)'],p.x,'linear',0)';
                else
                    lin2 = interp1([lend2,group2.x_nodes,1],[0,umv(group2.i_nodes)',0],p.x,'linear',0)';
                end
                lin = lin1 + lin2;
                
                lin(end) = [];
                
            else
                %No gap in the car
                if(v.i_nodes(1) ==0)
                    lin = interp1(v.x_nodes(2:end),umv(v.i_nodes(2:end)),p.x,'linear',0)';
                else
                    lin = interp1(v.x_nodes,umv(v.i_nodes),p.x,'linear',0)';
                    
                    
                end
                lin(end) = [];
            end
        otherwise
            lin = zeros(p.N-1,1);
            
    end
else
    if(v.x_nodes~=0)
        switch v.type
            case 1
                lin = interp1([0,v.x_nodes,1], [0,umv(v.i_nodes)',0], p.x, 'linear',0)';
                lin(end) = [];
            case {-1,2,4}
                lin = zeros(p.N,1);
                lin(v.i_nodes) = umv(v.i_nodes);
                lin(end) = [];
            case 3
                lend = v.i_nodes - ceil(v.length_car/2);
                rend = v.i_nodes + ceil(v.length_car/2);
                if(lend<1)
                    m = unique([0, p.x(lend+p.N), 1,1+p.x(rend),2]);
                    linlin = interp1([m, 1+v.x_nodes],[zeros(length(m),1)',umv(v.i_nodes)'], [p.x,p.x+1],'linear',0)';
                    lin = linlin(1:p.N) + linlin(p.N+1:end);
                elseif(rend>p.N+1)
                    m = unique([0, p.x(lend), 1,1+p.x(rend-p.N),2]);
                    linlin = interp1([m, v.x_nodes],[zeros(length(m),1)',umv(v.i_nodes)'], [p.x,p.x+1],'linear',0)';
                    lin = linlin(1:p.N) + linlin(p.N+1:end);
                elseif(lend==1&&rend~=p.N+1)
                    lin = interp1([0,v.x_nodes, p.x(rend), 1],[0,umv(v.i_nodes),0,0], [p.x],'linear',0)';
                elseif(rend==p.N+1&&lend~=1)
                    lin = interp1([0,p.x(lend),v.x_nodes, 1],[0,0,umv(v.i_nodes),0], [p.x],'linear',0)';
                elseif(rend==p.N+1&&lend==1)
                    lin = interp1([0,v.x_nodes, 1],[0,umv(v.i_nodes),0],[p.x],'linear',0)';
                else
                    lin = interp1([0, p.x(lend),v.x_nodes, p.x(rend), 1],[0, 0,umv(v.i_nodes),0,0], [p.x],'linear',0)';
                    
                end
                lin(end) = [];
        end
%         lin(end) = [];
    else
        lin = 0.*umv;
        lin(end) = [];
    end
    % lin = 1;
end
end

function [var] = Update_Overhaul(p, v)
% var = v;
% v.x_nodes = unique(floor(mod(v.x_nodes + v.velocity*p.t_adj,1)./p.dx).*p.dx);
switch v.type
    case -1
        v.i_nodes = unique(mod(v.i_nodes+v.direction*v.velocity,p.N));
%         if(v.i_nodes(1)==0)
%             v.i_nodes(1) = [];
%             v.i_nodes = [v.i_nodes, p.N];
%         end
        v.x_nodes = p.x(v.i_nodes);
    case {2, 4}
%         v.x_nodes = unique(floor(mod(v.x_nodes + v.direction*v.velocity*p.dt,1)./p.dx).*p.dx);
        v.x_nodes = unique(floor(mod(v.x_nodes_o + v.direction*v.velocity*p.t_adj,1)./p.dx).*p.dx);

        [maximum,index] = max(v.x_nodes(2:end)-v.x_nodes(1:end-1));
        if(maximum > p.dx*v.int_nodes)
%             v.i_nodes(1:index) = v.i_nodes(1:index);
            v.i_nodes = (v.x_nodes+p.dx)./p.dx;

        elseif(maximum == p.dx*2) %Gap caused when crossing over boundary
            v.x_nodes(1:index) = v.x_nodes(1:index) - p.dx;
            v.i_nodes = (v.x_nodes+p.dx)./p.dx;
            
        else
            v.i_nodes = (v.x_nodes+p.dx)./p.dx;
        end
        v.i_nodes = sort(v.i_nodes);
        v.x_nodes = sort(v.x_nodes);
%         if(v.i_nodes(1)==0)
%             v.i_nodes(1) = [];
%             v.i_nodes = [v.i_nodes];
%             v.x_nodes(1) = [];
%             v.x_nodes = [v.x_nodes, 1];
%         end
    case 3
        v.x_nodes = unique(floor(mod(v.x_nodes + v.direction*v.velocity*p.dt,1)./p.dx).*p.dx);
        v.i_nodes = (v.x_nodes+p.dx)./p.dx;
        if(v.x_nodes>=p.N+1)
            v.x_nodes = [0];
            v.i_nodes = [1];
        end
%         if(v.i_nodes(1)==0)
%             v.i_nodes(1) = [];
%             v.i_nodes = [v.i_nodes, p.N];
%             v.x_nodes(1) = [];
%             v.x_nodes = [v.x_nodes, 1];
%         end
    otherwise
        %         var = v;
end
var = v;
end
function [check] = Graphing_Overhaul(p,u,v)
persistent grapher_var fig grapher_u saver newGrapher
persistent spectrum
if(p.t == 0)
    grapher_var = [];
    fig = [];
    grapher_u = [];
    saver = [];
    newGrapher = [];
    spectrum = [];

end
m = length(v);
if(~isempty(fig)&&~isvalid(fig)&&isempty(grapher_u))
    clear fig;
end
if(isempty(fig))%||~isvalid(fig))
%     set(gca,'fontsize', 16);
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
    fig = gcf;
%     newGrapher = 
    
%     fig = figure('DefaultAxesFontSize',16);
% Set up figure properties:
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
%     set(gcf, 'Position', get(0, 'Screensize'));
%         figure('units','normalized','outerposition',[0 0 1 1])




else
    if(~ishghandle(fig))
        clear grapher_var fig grapher_u

        check  = false;
        return
    end 
end
if(gcf~=fig)
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)

    figure(fig)
    subplot(1,2,1);
end
subplot(1,2,1);
if(isempty(grapher_u))
    subplot(1,2,1);
    grapher_u.graph = plot(p.x,u.u, 'LineWidth',1);
    set(grapher_u.graph, 'YDataSource','u.u');
    hold on;

else
    subplot(1,2,1);
    refreshdata(grapher_u.graph,'caller');
end
if(isempty(grapher_var))
    subplot(1,2,1);
    grapher_var = repelem(struct,m,1);
    legendInfo = cell(1,2*m + 1);
    legendInfo{1} = 'Reference Solution';
    for i = 1:m
        switch v(i).type
            case 0
                grapher_var(i).graph = plot(p.x,v(i).u);
                set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).u',i));
                grapher_var(i).DAgraph = [];
                hold on;
                legendInfo{2*(i-1)+2} = 'Reference Solution';
                legendInfo{2*(i-1)+3} = '';
            case 1
                if(i==1)
                    grapher_var(i).graph = plot(p.x,v(i).v,'-g','LineWidth',1);
                    set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).v',i));
                    grapher_var(i).DAgraph = scatter(v(i).x_nodes,zeros(1,v(i).int_nodes),'ob');
                    set(grapher_var(i).DAgraph, 'XDataSource',sprintf('v(%i).x_nodes',i));
                    legendInfo{2*(i-1)+2} = 'Uniform Grid Solution';
                    legendInfo{2*(i-1)+3} = 'Uniform Gridpoints';
                    
                else
                    grapher_var(i).graph = plot(p.x,v(i).v,'--r','LineWidth',1);
                    set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).v',i));
                    grapher_var(i).DAgraph = scatter(v(i).x_nodes,zeros(1,v(i).int_nodes),'+r');
                    set(grapher_var(i).DAgraph, 'XDataSource',sprintf('v(%i).x_nodes',i));
                    legendInfo{2*(i-1)+2} = 'Uniform Grid Solution';
                    legendInfo{2*(i-1)+3} = 'Uniform Gridpoints';
                end

               
                
                hold on;
            case {-1, 2, 3, 4}
                grapher_var(i).graph = plot(p.x,v(i).v);
                set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).v',i));
                
                
                grapher_var(i).DAgraph = scatter(v(i).x_nodes,zeros(1,v(i).int_nodes),'+r');
                set(grapher_var(i).DAgraph, 'XDataSource',sprintf('v(%i).x_nodes',i));
                
                legendInfo{2*(i-1)+2} = 'Car Solution';
                legendInfo{2*(i-1)+3} = 'Car Gridpoints';
                hold on;
            otherwise
                
        end
    end
    y_amp = max([1.01/sqrt(p.beta/p.alpha), 1e-2]);
    axis tight manual;
    axis([0 p.L -(y_amp*1.2) y_amp*1.2]);
    title ( sprintf ('u(x, %1.3f)',p.t ));
    xlabel('X');
    ylabel('U(X)');
%     legend(legendInfo);
%     if(p.t_adj == 0.110 || p.t >= 3.9 && p.t <4.0 || p.t >= 8.3 && p.t < 8.4)
%         p.t
%     end
%     if(m~=0)
%         legend('Reference Solution','Uniform Solution','Uniform Gridpoints','Probe Solution','Probe Gridpoints');
%     else
%          legend('Reference Solution');
%     end
else
    for i = 1:m
        refreshdata(grapher_var(i).graph,'caller');
        if(~isempty(grapher_var(i).DAgraph))
            refreshdata(grapher_var(i).DAgraph,'caller');
        end
    end
    title ( sprintf ('u(x, %1.3f)',p.t_adj ));
%     if(p.t_adj >= 1.27)
%         pause
%     end
    
end
subplot(1,2,2);
u_hat = fft(u.u);
loglog(1:p.N/2, abs(u_hat(1:p.N/2)));
hold on;
loglog([p.N/4,p.N/4],[1e-20,1e4])
hold off;
axis([0 p.N/2 1e-18 1e4]);

drawnow; 
% savecheck = true;
% filename = 'carexample2.gif';
% % Capture the plot as an image
% frame = getframe(fig);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% 
% % Write to the GIF File
% if  isempty(saver)&&savecheck
%     imwrite(imind,cm,filename,'gif','DelayTime',0.1,  'Loopcount',inf);
%     saver = true;
% elseif(savecheck)
%     imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'WriteMode','append');
% end
check = true;
return
end