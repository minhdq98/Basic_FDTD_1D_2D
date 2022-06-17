%% 1D FDTD simulation program for dispersive materials with PML layers
% Written by Quang Minh Dinh

clear 

%% I. Set up the simulation

% 1. Basic parameters
L = 80; % Length of the simulation space, unit = um
Nx = 800; % Number of grid points
dx = L/Nx; % Size of a pixel (distance between two adjacent grid points)
c0 = 300; % Free-space speed of light, unit = um/ps
dt = 0.5/(c0*sqrt(1/dx^2)); % Time step that meets the CFL condition
eps0 = 8.85e-6; % Epsilon_0, unit: F/um
mu0 = 1.26; % Mu_0, unit: H/um
x = linspace(0,L,Nx); % Coordinates of the grid points
omega = 2*pi*100; % Angular frequency of waves that we are investigating, unit: 1/ps

% 2. Customize the environment
mu = ones(1,Nx); % Permeability mu = 1 everywhere
eps = ones(1,Nx); % Preset the permittivity eps = 1 everywhere initially
sigma = zeros(1,Nx); % Preset the conductivity sigma = 0 everywhere intially
% Depending on the chosen scenario, the permittivity and conductivity
% will be given accordingly
scen = 1; % Selected scenario
if scen == 1
    % Insert a near-zero-index slab of thickness N1 at x1
    N1 = Nx/4;
    x1 = Nx/2;
    % Material properties
    eps_inf = 1; % Permittivity at very high frequency
    omegap = omega*1.01; % Plasma frequency, unit: 1/ps
    omega0 = omega*0.1; % Resonant frequency, unit: 1/ps
    gamma = 0; % Scattering rate, unit: 1/ps
    mu_inf = 1; % Permeability at very high frequency
    omegapm = omega*1.01; % Magnetic plasma frequency, unit: 1/ps
    omega0m = omega*0.1; % Magnetic resonant frequency, unit: 1/ps
    gammam = 0; % Magnetic scattering rate, unit: 1/ps
end
if scen == 2
    % Insert a double-negative slab of thickness N2 at x2
    N1 = Nx/4;
    x1 = Nx/2;
    % Material properties
    eps_inf = 1; % Permittivity at very high frequency
    omegap = omega*1.2; % Plasma frequency, unit: 1/ps
    omega0 = omega*0.3; % Resonant frequency, unit: 1/ps
    gamma = 0; % Scattering rate, unit: 1/ps
    mu_inf = 1; % Permeability at very high frequency
    omegapm = omega*1.2; % Magnetic plasma frequency, unit: 1/ps
    omega0m = omega*0.3; % Magnetic resonant frequency, unit: 1/ps
    gammam = 0; % Magnetic scattering rate, unit: 1/ps
end
if scen == 3
    % Insert a normal dispersive dielectric slab of thickness N2 at x2
    N1 = Nx/4;
    x1 = Nx/2;
    % Material properties
    eps_inf = 1; % Permittivity at very high frequency
    omegap = omega*5; % Plasma frequency, unit: 1/ps
    omega0 = omega*4; % Resonant frequency, unit: 1/ps
    gamma = 0; % Scattering rate, unit: 1/ps
    mu_inf = 1; % Permeability at very high frequency
    omegapm = 0; % Magnetic plasma frequency, unit: 1/ps
    omega0m = 0; % Magnetic resonant frequency, unit: 1/ps
    gammam = 0; % Magnetic scattering rate, unit: 1/ps
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
M = -dt./(mu0*dx*(mu+1/2*dt*sigma_star));
C = eps*eps0./(eps*eps0+sigma*dt);
N = -(dt/dx)./(eps*eps0+sigma*dt);

% 5. Initialize the fields to be zero everywhere
E = zeros(1,Nx);
H = zeros(1,Nx);
D = zeros(1,Nx); % Displacement field
B = zeros(1,Nx); % Magnetic flux density field
P = zeros(1,Nx); % Polarization density field
J = zeros(1,Nx); % Magnetic polarization density field
P1 = zeros(1,Nx); % Polarization density at one previous iteration
P2 = zeros(1,Nx); % Polarization density at two previous iteration
J1 = zeros(1,Nx); % Magnetic polarization density at one previous iteration
J2 = zeros(1,Nx); % Magnetic polarization density at two previous iteration

%% II. Run simulation and display results

% 1. Basic parameters
i = 0; % Time step, count the number of iterations
n = 2; % Ignore 2 grid points adjacent to each boundary
% The region of the dispersive material
NF = x1-N1/2;
NB = x1+N1/2;
% The region which is simulated
N3 = n;
N4 = Nx-n;

% 2. Source parameters
Ns = NPML+1; % Position of the source
T = 20; % Settling time of the source
A0 = 1; % Stable amplitude of the source

% 3. Start the simulation
while(1) % Keep the simulation running until the user cancels it
    % 3.1. Insert the source, which is sinusoidal but initially increases 
    % from 0 and gradually reaches a stable amplitude after a settling time
    E(1,Ns) = E(1,Ns)+A0*(1-exp(-(i/T)^2))*cos(omega*i*dt);
    % 3.2. Update the fields in free-space and in the PML
    % In front of the dispersive medium
    H(1,N3:NF-1) = A(1,N3:NF-1).*H(1,N3:NF-1)+M(1,N3:NF-1).*(E(1,N3+1:NF)-E(1,N3:NF-1));
    E(1,N3:NF-1) = C(1,N3:NF-1).*E(1,N3:NF-1)+N(1,N3:NF-1).*(H(1,N3:NF-1)-H(1,N3-1:NF-2));
    % Behind the dispersive medium
    H(1,NB+1:N4) = A(1,NB+1:N4).*H(1,NB+1:N4)+M(1,NB+1:N4).*(E(1,NB+2:N4+1)-E(1,NB+1:N4));
    E(1,NB+1:N4) = C(1,NB+1:N4).*E(1,NB+1:N4)+N(1,NB+1:N4).*(H(1,NB+1:N4)-H(1,NB:N4-1));
    % 3.3. Update the B, J, and H fields in the dispersive medium
    B(1,NF:NB) = B(1,NF:NB)-(dt/dx)*(E(1,NF+1:NB+1)-E(1,NF:NB));
    J(1,NF:NB) = (2-dt^2*omega0m^2)*J1(1,NF:NB)/(1+0.5*dt*gammam)+...
        (0.5*dt*gamma-1)*J2(1,NF:NB)/(1+0.5*dt*gammam)+...
        dt^2*mu0*omegapm^2*H(1,NF:NB)/(1+0.5*dt*gammam);
    H(1,NF:NB) = (B(1,NF:NB)-J(1,NF:NB))/(mu0*mu_inf);
    % 3.4. Update the D, P, and E fields in the dispersive medium
    D(1,NF:NB) = D(1,NF:NB)-(dt/dx)*(H(1,NF:NB)-H(1,NF-1:NB-1));
    P(1,NF:NB) = (2-dt^2*omega0^2)*P1(1,NF:NB)/(1+0.5*dt*gamma)+...
        (0.5*dt*gamma-1)*P2(1,NF:NB)/(1+0.5*dt*gamma)+...
        dt^2*eps0*omegap^2*E(1,NF:NB)/(1+0.5*dt*gamma);
    E(1,NF:NB) = (D(1,NF:NB)-P(1,NF:NB))/(eps0*eps_inf);
    % 3.5. Store the P and J fields to use in the next iterations
    P2(1,NF:NB) = P1(1,NF:NB);
    P1(1,NF:NB) = P(1,NF:NB);
    J2(1,NF:NB) = J1(1,NF:NB);
    J1(1,NF:NB) = J(1,NF:NB);
    % 3.7. Display the results
    figure(1)
    plot(x,real(E),'b');
    axis([0 L -2 2]);
    i = i+1;
end
