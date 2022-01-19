% function [f1_hat, f2_hat, g_hat] = OlsonTitiForcing(Nx,Ny,G,nu)
function g_hat = OlsonTitiForcing(Nx,Ny,nu,G)
% Author: Adam Larios. Last edit: 2020.12.17
% Calculate the function that you put on the stream function formulation,
% before inverse Laplacian is applied.  That is, g is at vorticity level.
% Reference: Eric Olson, Edriss S. Titi. "Determining Modes and Grashof Number in
%            2D Turbulence", Theoretical and Computational Fluid Dynamics, Vol. 22, No. 5, 2008, pp. 327â339.

% From the paper:
% Fourier modes for g_0 = nabla Ã f_0. Since the reality
% condition implies Ä?_{âk} = Ä?_k we list only half the modes here.

% Therefore, need symmetrize, and then "uncurl" g.
% To uncurl, we set
% -Laplacian f1 = -dg/dy
% -Laplacian f2 =  dg/dx
% and solve for f1 and f2.

% G = desired Grashof number for viscosity nu;
% G = ||f||_{L^2}*(L/(2*pi*nu))^2,
% where L = 2*pi = side length of the periodic box.


if ~(exist('Nx','var'))
    Nx = 2^9;  % Number of x gridpoints (uniform).
end
if ~(exist('Ny','var'))
    Ny = Nx; % Number of y gridpoints (uniform).
end
if (Nx < 32) || (Ny < 32)
    error('Need Nx>=32 and Ny>=32 for Olson-Titi forcing');
end
if ~(exist('G','var'))
%     G = 250000;  % Target Grashof number
%     G = 1e6;
%     G = 1e6; % Gesho-Olson-Tit paper
G = 2.5e6;
end
if ~(exist('nu','var'))
    nu = 0.0001;  % viscosity
end

Lx = 2*pi; % Length of physical domain in x direction.
Ly = 2*pi; % Length of physical domain in y direction.
parseval = sqrt(Lx*Ly)/Nx/Ny;
lambda1 = min([(2*pi/Lx)^2,(2*pi/Ly)^2]); % Eigenfunctions sin(2*pi*x/Lx) and so on.

% f_hat = zeros(Nx,Ny);
% 
% f_hat(5,12) = -1.68185 + -1.37884*1i;
% f_hat(4,12) = -0.501504 + -4.39707*1i;
% f_hat(Nx-10,2) = 3.22868 + 3.57232*1i;
% f_hat(11,2) = -0.168267 + 3.63812*1i;
% f_hat(Nx-10,3) = 1.83238 + 4.56593*1i;
% f_hat(11,3) = 1.08206 + -2.65437*1i;
% f_hat(Nx-10,4) = 1.23101 + 2.32092*1i;
% f_hat(11,4) = -4.36393 + 4.81623*1i;
% f_hat(Nx-10,5) = 0.772188 + -2.22768*1i;
% f_hat(11,5) = 0.352859 + 3.8498*1i;
% f_hat(Nx-9,6) = 2.95066 + -3.88721*1i;
% f_hat(10,6) = -4.65628 + -0.853542*1i;
% f_hat(Nx-9,7) = 2.14982 + 3.9452*1i;
% f_hat(Nx-7,7) = 0.847957 + 6.09509*1i;
% f_hat(10,7) = -0.933718 + 1.34388*1i;
% f_hat(Nx-8,8) = 3.75659 + -1.43354*1i;
% f_hat(9,8) = -3.07828 + 1.434*1i;
% f_hat(Nx-7,9) = 1.57689 + -3.96814*1i;
% f_hat(Nx-5,9) = 0.689446 + -2.66231*1i;
% f_hat(8,9) = 2.87397 + -1.28323*1i;
% f_hat(Nx-6,10) = 3.06962 + 1.30303*1i;
% f_hat(Nx-4,10) = 1.3582 + -2.86109*1i;
% f_hat(7,10) = -1.08488 + -0.661055*1i;
% f_hat(Nx-5,11) = 0.178138 + -2.35925*1i;
% f_hat(Nx-3,11) = -1.30734 + -1.20256*1i;
% f_hat(Nx-1,11) = -4.54014 + 2.54287*1i;
% f_hat(1,11) = 2.28551 + -4.28632*1i;
% f_hat(3,11) = -2.75146 + 0.206303*1i;
% f_hat(5,11) = 0.0705579 + -0.807738*1i;
% f_hat(7,11) = 1.97628 + -1.26563*1i;
% f_hat(Nx-2,12) = 4.63344 + 3.9968*1i;
% f_hat(1,12) = 3.03338 + -1.77921*1i;
% f_hat(2,12) = -3.64792 + 6.35805*1i;
% f_hat(11,1) = 1.24568 + 0.999582*1i;
% f_hat(12,1) = -2.45397 + -2.1001*1i;
% f_hat(Nx-9,2) = -3.6533 + 4.34302*1i;
% f_hat(12,2) = 1.74878 + 0.537828*1i;
% f_hat(Nx-9,3) = -2.27939 + 5.32842*1i;
% f_hat(12,3) = -0.12495 + -0.439438*1i;
% f_hat(Nx-9,4) = 3.24541 + -2.56086*1i;
% f_hat(12,4) = -3.25672 + -2.04313*1i;
% f_hat(Nx-9,5) = -1.18464 + 3.6144*1i;
% f_hat(12,5) = -0.0779563 + -0.483413*1i;
% f_hat(Nx-8,6) = 1.43467 + -3.13185*1i;
% f_hat(11,6) = -0.678476 + 1.56749*1i;
% f_hat(Nx-8,7) = 0.794682 + -0.405619*1i;
% f_hat(9,7) = 0.412921 + 1.72437*1i;
% f_hat(11,7) = -0.960698 + 1.73383*1i;
% f_hat(Nx-7,8) = -2.33465 + -1.91763*1i;
% f_hat(10,8) = -1.4385 + 2.08662*1i;
% f_hat(Nx-6,9) = -0.206682 + 2.31327*1i;
% f_hat(7,9) = -0.956604 + -2.33234*1i;
% f_hat(9,9) = 0.140247 + 0.118369*1i;
% f_hat(Nx-5,10) = 2.42145 + -0.448085*1i;
% f_hat(6,10) = -3.58454 + 1.92307*1i;
% f_hat(8,10) = 2.62845 + 3.4255*1i;
% f_hat(Nx-4,11) = 1.83893 + 3.99093*1i;
% f_hat(Nx-2,11) = 3.46985 + 0.165801*1i;
% f_hat(1,11) = -1.81676 + -6.38674*1i;
% f_hat(2,11) = 4.89391 + -0.112915*1i;
% f_hat(4,11) = -2.18555 + -0.921325*1i;
% f_hat(6,11) = 1.54712 + -1.76146*1i;
% f_hat(Nx-3,12) = -4.18572 + -0.412749*1i;
% f_hat(Nx-1,12) = 0.161582 + 4.15786*1i;
% f_hat(1,12) = -0.608278 + -1.09338*1i;
% f_hat(3,12) = 1.35454 + -1.92183*1i;

% Data directly from the paper.
k1 = [4,3,-11,10,-11,10,-11,10,-11,10,-10,9,-10,-8,9,-9,8,-8,-6,7,-7,-5,6,-6,-4,-2,0,2,4,6,-3,-1,1,10,11,-10,11,-10,11,-10,11,-10,11,-9,10,-9,8,10,-8,9,-7,6,8,-6,5,7,-5,-3,-1,1,3,5,-4,-2,0,2];
k2 = [11,11,1,1,2,2,3,3,4,4,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,10,10,10,10,10,11,11,11,0,0,1,1,2,2,3,3,4,4,5,5,6,6,6,7,7,8,8,8,9,9,9,10,10,10,10,10,10,11,11,11,11];
gk1 = 1e-4*[-1.68185,-0.501504,3.22868,-0.168267,1.83238,1.08206,1.23101,-4.36393,0.772188,0.352859,2.95066,-4.65628,2.14982,0.847957,-0.933718,3.75659,-3.07828,1.57689,0.689446,2.87397,3.06962,1.3582,-1.08488,0.178138,-1.30734,-4.54014,2.28551,-2.75146,0.0705579,1.97628,4.63344,3.03338,-3.64792,1.24568,-2.45397,-3.6533,1.74878,-2.27939,-0.12495,3.24541,-3.25672,-1.18464,-0.0779563,1.43467,-0.678476,0.794682,0.412921,-0.960698,-2.33465,-1.4385,-0.206682,-0.956604,0.140247,2.42145,-3.58454,2.62845,1.83893,3.46985,-1.81676,4.89391,-2.18555,1.54712,-4.18572,0.161582,-0.608278,1.35454];
gk2 = 1e-4*[-1.37884,-4.39707,3.57232,3.63812,4.56593,-2.65437,2.32092,4.81623,-2.22768,3.8498,-3.88721,-0.853542,3.9452,6.09509,1.34388,-1.43354,1.434,-3.96814,-2.66231,-1.28323,1.30303,-2.86109,-0.661055,-2.35925,-1.20256,2.54287,-4.28632,0.206303,-0.807738,-1.26563,3.9968,-1.77921,6.35805,0.999582,-2.1001,4.34302,0.537828,5.32842,-0.439438,-2.56086,-2.04313,3.6144,-0.483413,-3.13185,1.56749,-0.405619,1.72437,1.73383,-1.91763,2.08662,2.31327,-2.33234,0.118369,-0.448085,1.92307,3.4255,3.99093,0.165801,-6.38674,-0.112915,-0.921325,-1.76146,-0.412749,4.15786,-1.09338,-1.92183];

gk = gk1 + 1i*gk2;

% Adjust to Matlab-style indices
for  i = 1:length(k1)
    if k1(i) < 0
        k1(i) = Nx + k1(i) +1;
    else
        k1(i) = k1(i)  + 1;
    end
end
k2 = k2 + 1;  % Note that all k2's are non-negative here.

g_hat = zeros(Nx,Ny);
for i = 1:length(k1)    
    if (k1(i) > 1) && (k2(i) > 1)
        g_hat(k1(i),k2(i)) = gk(i);
        g_hat(Nx - k1(i) +2  ,Ny - k2(i) + 2) = conj(gk(i));
    end
    
    if k1(i) == 1 % We are on the real axis.
        g_hat(1,k2(i)) = gk(i);
        g_hat(1,Ny - k2(i) + 2) = conj(gk(i));
    end
    if k2(i) == 1 % We are on the real axis.
        g_hat(k1(i),1) = gk(i);
        g_hat(Ny - k1(i) + 2,1) = conj(gk(i));
    end
    
end

% g = ifft2(g_hat);
% subplot(1,3,1)
% pcolor(abs(fftshift(g_hat))); axis('square');
% subplot(1,3,2)
% pcolor(abs((g_hat))); axis('square');
% subplot(1,3,3)
% pcolor(g); axis('square');
% shading interp

%% === Next, we uncurl to obtain f ========================================

% Wave numbers k (we also multiply k by i for convenience).
ikx    = (1i*[0:Nx/2-1, 0, -Nx/2+1:-1]*(2*pi/Lx)).'; % Note: .' is the non-conjugate transpose
kx_sq  =    ([0:Nx/2,      -Nx/2+1:-1].^2*(2*pi/Lx)^2).';

% Wave numbers k (we also multiply k by i for convenience).
iky   = 1i*[0:Ny/2-1, 0, -Ny/2+1:-1]*(2*pi/Ly);
ky_sq =    [0:Ny/2,      -Ny/2+1:-1].^2*(2*pi/Ly)^2;

% Wave numbers of Laplacian operator.
k_lap_inv = zeros(Nx,Ny);
for j = 1:Ny
    for i = 1:Nx
        norm_k_sq = (kx_sq(i) + ky_sq(j));
%         k_lap = -(kx_sq(i) + ky_sq(j));
        k_lap_inv(i,j) = -1/norm_k_sq;
    end
end
k_lap_inv(abs(k_lap_inv) == Inf) = 0;
k_lap_inv(isnan(k_lap_inv)) = 0;
k_lap_inv(1) = 0;

% %%% Mike (Nov 10, 2020) says: 
% %fnrm=2\pi  sqrt(sum_k  |g_k|^2/|k|^2 ),
% %f=G*nu**2*f/fnrm
% fnrm = 0;
% for j = 1:Ny
%     for i = 1:Nx
%         norm_k_sq = (kx_sq(i) + ky_sq(j));
%         if norm_k_sq ~= 0
%             fnrm = fnrm + abs(g_hat(i,j))^2/norm_k_sq;
%         end
%     end
% end
% fnrm=2*pi*sqrt(fnrm); % Should be 0.0025 as in Olson-Titi-2008
% 
% % % ==============================================
% % g = curl(f); curl(g) = -Laplace(f); f = -invLap*curl(g)
% % Also:
% % psi_f = forcing at stream function level; that is, f = curl(psi_f)
% % So:
% % g = curl(f) = -Laplace(psi_f); so psi_f = -invLap*g
% 
% % % To compute all of this, we follow these steps:
% % f1_hat = -(G*nu^2/fnrm)*k_lap_inv.*bsxfun(@times,p.iky,g_hat);
% % f2_hat =  (G*nu^2/fnrm)*k_lap_inv.*bsxfun(@times,p.ikx,g_hat);
% % g_hat = bsxfun(@times,p.ikx,f2_hat) - bsxfun(@times,p.iky,f1_hat);
% % % However, this can be simplified, since curl(curl) = -Laplace
% g_hat = (G*nu^2/fnrm)*g_hat;
% % % If we need to get back to stream-function forcing, we can do this:
% % psi_f_hat = -k_lap_inv.*g_hat;
% % % ==============================================

% norm(g_hat)*sqrt(dx*dy)

% % Compute f from g as follows:
%     g := curl(f)
%     curl(g) = curl(curl(f)) = -Lap(f)
%     f = -inv_Lap(curl(g))
f1_hat = (-k_lap_inv).*(-bsxfun(@times,iky,g_hat));
f2_hat = (-k_lap_inv).*( bsxfun(@times,ikx,g_hat));
f_L2norm = sqrt(norm(f1_hat,'fro')^2 + norm(f2_hat,'fro')^2)*parseval;

% fprintf('G = %g\n',norm(f,'fro')*sqrt((2*pi/512)*(2*pi/512))/nu^2);

current_G = f_L2norm/(lambda1*nu^2);
% Rescale to get target Grashof number.
f1_hat = (G/current_G)*f1_hat;
f2_hat = (G/current_G)*f2_hat;

g_hat = bsxfun(@times,ikx,f2_hat) - bsxfun(@times,iky,f1_hat);

% norm(g_hat,'fro')*parseval

% f_L2norm = sqrt(norm(f1_hat,'fro')^2 + norm(f2_hat,'fro')^2)*parseval;
% f_L2norm/(lambda1*nu^2)
% 
% 555
% g_hat = (G/current_G)*g_hat;




% g_x_hat =  bsxfun(@times,ikx,g_hat);
% g_y_hat =  bsxfun(@times,iky,g_hat);
% 
% f1_hat = -g_y_hat;
% f2_hat =  g_x_hat;
% 
% % f1 = ifft2(f1_hat);
% % f2 = ifft2(f2_hat);
% 

% norm = (max([max(abs(f1_hat(:))),max(abs(f2_hat(:)))]));
% f_L2norm = sqrt(norm(f1_hat(:))^2 + norm(f2_hat(:))^2)*parseval;
% g_L2norm = norm(g_hat,'fro')*parseval;
% % f_L2norm/g_L2norm
% 
% f1_hat = f1_hat*(G*nu^2/f_L2norm);
% f2_hat = f2_hat*(G*nu^2/f_L2norm);

end

%% List from original paper in plain text.
%   Uncomment and pipe this to awk to get the above list 
%   Note the shift in index, since Matlab starts at 1.
%   We also pipe to sed to account for the negative and zeroindices.
%   xsel | awk '{print "f_hat("$1+1","$2+1") = "$3" + "$4"*1i;"}' | sed 's/(-/(Nx-/' | sed 's/(0/(1/'

% 4 11 -1.68185 -1.37884
% 3 11 -0.501504 -4.39707
% -11 1 3.22868 3.57232
% 10 1 -0.168267 3.63812
% -11 2 1.83238 4.56593
% 10 2 1.08206 -2.65437
% -11 3 1.23101 2.32092
% 10 3 -4.36393 4.81623
% -11 4 0.772188 -2.22768
% 10 4 0.352859 3.8498
% -10 5 2.95066 -3.88721
% 9 5 -4.65628 -0.853542
% -10 6 2.14982 3.9452
% -8 6 0.847957 6.09509
% 9 6 -0.933718 1.34388
% -9 7 3.75659 -1.43354
% 8 7 -3.07828 1.434
% -8 8 1.57689 -3.96814
% -6 8 0.689446 -2.66231
% 7 8 2.87397 -1.28323
% -7 9 3.06962 1.30303
% -5 9 1.3582 -2.86109
% 6 9 -1.08488 -0.661055
% -6 10 0.178138 -2.35925
% -4 10 -1.30734 -1.20256
% -2 10 -4.54014 2.54287
% 0 10 2.28551 -4.28632
% 2 10 -2.75146 0.206303
% 4 10 0.0705579 -0.807738
% 6 10 1.97628 -1.26563
% -3 11 4.63344 3.9968
% -1 11 3.03338 -1.77921
% 1 11 -3.64792 6.35805
% 10 0 1.24568 0.999582
% 11 0 -2.45397 -2.1001
% -10 1 -3.6533 4.34302
% 11 1 1.74878 0.537828
% -10 2 -2.27939 5.32842
% 11 2 -0.12495 -0.439438
% -10 3 3.24541 -2.56086
% 11 3 -3.25672 -2.04313
% -10 4 -1.18464 3.6144
% 11 4 -0.0779563 -0.483413
% -9 5 1.43467 -3.13185
% 10 5 -0.678476 1.56749
% -9 6 0.794682 -0.405619
% 8 6 0.412921 1.72437
% 10 6 -0.960698 1.73383
% -8 7 -2.33465 -1.91763
% 9 7 -1.4385 2.08662
% -7 8 -0.206682 2.31327
% 6 8 -0.956604 -2.33234
% 8 8 0.140247 0.118369
% -6 9 2.42145 -0.448085
% 5 9 -3.58454 1.92307
% 7 9 2.62845 3.4255
% -5 10 1.83893 3.99093
% -3 10 3.46985 0.165801
% -1 10 -1.81676 -6.38674
% 1 10 4.89391 -0.112915
% 3 10 -2.18555 -0.921325
% 5 10 1.54712 -1.76146
% -4 11 -4.18572 -0.412749
% -2 11 0.161582 4.15786
% 0 11 -0.608278 -1.09338
% 2 11 1.35454 -1.92183