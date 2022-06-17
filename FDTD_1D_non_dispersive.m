%% 1D FDTD simulation program for non-dispersive, non-magnetic materials with PML layers
% Written by Quang Minh Dinh

clear

%% I. Set up the simulation

% 1. Basic parameters
L = 80; % Length of the simulation space, unit = um
Nx = 800; % Number of grid points
dx = L/Nx; % Size of a pixel (distance between two adjacent grid points)
c0 = 300; % Free-space speed of light, unit = um/ps
dt = 1/(c0*sqrt(1/dx^2)); % Time step that meets the CFL condition
eps0 = 8.85e-6; % Epsilon_0, unit: F/um
mu0 = 1.26; % Mu_0, unit: H/um
x = linspace(0,L,Nx); % Coordinates of the grid points

% 2. Customize the permittivity, permeability and conductivity
mu = ones(1,Nx); % Permeability mu = 1 everywhere
eps = ones(1,Nx); % Preset the permittivity eps = 1 everywhere initially
sigma = zeros(1,Nx); % Preset the conductivity sigma = 0 everywhere intially 
% Depending on the chosen scenario, the permittivity and conductivity
% will be given accordingly
scen = 4; % Chosen scenario
if scen == 1
    N1 = Nx/4;
    eps1 = 4;
    x1 = Nx/2;
    % Insert a lossless dielectric slab with permittivity eps1 and
    % thickness N1 pixels at x1
    eps(1,round(x1-N1/2):round(x1+N1/2)) = eps1; 
end
if scen == 2
    N2 = Nx/4;
    eps2 = 2;
    sigma2 = 0.2;
    x2 = Nx/2;
    % Insert a lossy dielectric slab with permittivity eps2, conductivity
    % sigma2 and thickness N2 pixels at x2
    eps(1,round(x2-N2/2):round(x2+N2/2)) = eps2; 
    sigma(1,round(x2-N2/2):round(x2+N2/2)) = sigma2;
end
if scen == 3
    N3 = Nx/10;
    sigma3 = 100;
    x3 = Nx/2;
    % Insert a conductor slab with conductivity sigma3 and thickness N3 at
    % x3
    sigma(1,round(x3-N3/2):round(x3+N3/2)) = sigma3;
end
if scen == 4
    N41 = Nx/200;
    N42 = Nx/100;
    eps41 = 4;
    eps42 = 1;
    x4 = Nx/2.5;
    % Insert a Bragg reflector with permittivity eps41, eps42, thickness
    % N41, N42, starting at x4
    eps(1,x4:x4+N41) = eps41;
    eps(1,x4+N41:x4+N41+N42) = eps42;
    eps(1,x4+N41+N42:x4+N41+N42+N41) = eps41;
    eps(1,x4+N41+N42+N41:x4+N41+N42+N41+N42) = eps42;
    eps(1,x4+N41+N42+N41+N42:x4+N41+N42+N41+N42+N41) = eps41;
    eps(1,x4+N41+N42+N41+N42+N41:x4+N41+N42+N41+N42+N41+N42) = eps42;
    eps(1,x4+N41+N42+N41+N42+N41+N42:x4+N41+N42+N41+N42+N41+N42+N41) = eps41;
end

% 3. Create the perfectly matched layers
NPML = Nx/10; % Thickness of the PMLs
r_required = 1e-40; % Required reflection coefficient
m = 3; % Polynomial order
sigma_max = -(m+1)*log(r_required)/(2*NPML*sqrt(mu0/eps0));
P = ((1:NPML)./NPML).^m*sigma_max; % Conductivity profile of the PML
% Put that conductivity profile at the right and left boundaries
sigma(1,Nx-NPML+1:Nx) = P; 
sigma(1,1:NPML) = fliplr(P);
% Magnetic conductivity
sigma_star = (sigma.*mu0)./(eps0);

% 4. Compute the coefficient matrices
A = (mu-1/2*dt*sigma_star)./(mu+1/2*dt*sigma_star);
B = -dt./(mu0*dx*(mu+1/2*dt*sigma_star));
C = eps*eps0./(eps*eps0+sigma*dt);
D = -(dt/dx)./(eps*eps0+sigma*dt);

% 5. Initialize the fields to be zero everywhere
E = zeros(1,Nx);
H = zeros(1,Nx);

%% II. Run simulation and display results

% 1. Basic parameters
i = 0; % Time step, count the number of iterations
n = 2; % Ignore 2 grid points adjacent to each boundary
N1 = n; % The fields are only updated in the region between N1, N2
N2 = Nx-n;

% 2. Source parameters
omega = 2*pi*95; % Angular frequency of the source, unit 1/ps
Ns = NPML+1; % Position of the source
T = 20; % Settling time of the source
A0 = 1; % Stable amplitude of the source

% 3. Start simulation
while(1) % Keep the simulation running until the user cancels it
    % Insert the source, which is sinusoidal but initially increases from 0
    % and gradually reaches a stable amplitude after a settling time
    E(1,Ns) = E(1,Ns)+A0*(1-exp(-(i/T)^2))*cos(omega*i*dt);
    % Updating the H and E fields
    H(1,N1:N2) = A(1,N1:N2).*H(1,N1:N2)+B(1,N1:N2).*(E(1,N1+1:N2+1)-E(1,N1:N2));
    E(1,N1:N2) = C(1,N1:N2).*E(1,N1:N2)+D(1,N1:N2).*(H(1,N1:N2)-H(1,N1-1:N2-1));
    % Display the results
    figure(1)
    yyaxis left
    plot(x,real(E),'b');
    axis([0 L -2 2]);
    yyaxis right
    ylim([0 2])
    plot(x,eps-ones(1,Nx));
    i = i+1;
end