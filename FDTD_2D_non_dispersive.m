%% 2D FDTD simulation program for non-dispersive, non-magnetic materials with PML layers
% Written by Quang Minh Dinh

clear

%% I. Set up the parameters

% 1. Basic parameters
L = 120; % Length of the simulation space, unit: um
W = 120; % Width of the simulation space, unit: um
Nx = 1200; % Number of grid points in the x axis
Ny = 1200; % Number of grid points in the y axis
dx = L/Nx; % Size of a pixel, unit: um
c0 = 300; % Speed of light in free-space, unit: um/ps
dt = 1/(c0*sqrt(2/dx^2)); % Time step that meets the CFL condition
eps0 = 8.85e-6; % Epsilon_0, unit: F/um
mu0 = 1.26; % Mu_0, unit: H/um
x = linspace(0,L,Nx); % List of coordinates in the x axis
y = linspace(0,W,Ny); % List of coordinates in the y axis

% 2. Create the permittivity, permeability, and conductivity profiles,
% which depend on the selected scenario
mu = ones(Ny,Nx); % Preset the permeability to be one everywhere
eps = ones(Ny,Nx); % Preset the permittivity to be one everywhere
sigma = zeros(Ny,Nx); % Preset the conductivity to be zero everywhere
scen = 4; % Chosen scenario
if scen == 1
    N1 = Ny/4;
    eps1 = 4;
    y1 = Ny/2;
    % Insert a dielectric slab of thickness N1, permittivity eps1 at y1
    eps(round(y1-N1/2):round(y1+N1/2),:) = eps1;
end
if scen == 2
    N2 = Ny/4;
    sigma2 = 0.0001;
    eps2 = 2;
    y2 = Ny/2;
    % Insert a lossy dielectric slab of thickness N2, permittivity eps2,
    % conductivity sigma2 at y2
    eps(round(y2-N2/2):round(y2+N2/2),:) = eps2;
    sigma(round(y2-N2/2):round(y2+N2/2),:) = sigma2;
end
if scen == 3
    N3 = Ny/10;
    sigma3 = 100;
    y3 = Ny/2;
    % Insert a conductor slab of thickness N3, conductivity sigma3 at y3
    sigma(round(y3-N3/2):round(y3+N3/2),:) = sigma3;
end
if scen == 4
    N4 = Ny/4;
    Lambda4 = Nx/20;
    eps4_var = 3;
    eps4_stat = 4;
    y4 = Ny/2;
    % Insert a Bragg gratings of Thickness N4, periodicity Gamma4 in the x
    % axis, permittivity varies between eps4_stat +/- eps4_var 
    % at y4
    eps(round(y4-N4/2):round(y4+N4/2),:) = repmat(eps4_stat+eps4_var*square(2*pi*x/(Lambda4*dx)),N4+1,1);
end

% 3. Create perfectly matched layers
NPML = Nx/10; % Thickness of the PMLs
r_required = 1e-40; % Required reflection coefficient
m = 3; % Polynomial order
sigma_max = -(m+1)*log(r_required)/(2*NPML*sqrt(mu0/eps0));
P = ((1:NPML)./NPML).^m*sigma_max; % Conductivity profile of the PML
% Put the conductivity profile at the left and right boundaries
sigma(:,Nx-NPML+1:Nx) = repmat(P,Ny,1); 
sigma(:,1:NPML) = fliplr(repmat(P,Ny,1));
% Put the conductivity profile at the top and bottom boundaries
sigma(1:NPML,:) = flip(repmat(P.',1,Nx));
sigma(Ny-NPML+1:Ny,:) = repmat(P.',1,Nx);
% Magnetic conductitivy
sigma_star = sigma*mu0/eps0;

% 4. Compute the coefficient matrices
A = (mu-1/2*dt*sigma_star)./(mu+1/2*dt*sigma_star);
B = dt./(mu0*dx*(mu+1/2*dt*sigma_star));
C = eps*eps0./(eps*eps0+sigma*dt);
D = dt./(dx*(eps*eps0+sigma*dt));

% We will have the Maxwell equations like this:
% Hz = A.*Hz-B.*(Ey-Ey2)+B.*(Ex-Ex2), Ex = C.*Ex+D.*(Hz2-Hz), Ey = C.*Ey+D.*(Hz-Hz2)

% 5. Initialize the fields to to zero everywhere
Hz = zeros(Ny,Nx);
Ex = zeros(Ny,Nx);
Ey = zeros(Ny,Nx);

%% II. Run simulation and display results

% 1. Basic parameters
i = 0; % Time step, count the number of iterations
n = 2; % Ignore n pixels adjacent to each boundaries
N1 = n; % The fields are only updated in the region between (N1,N3), (N1,N4), (N2,N4), (N2,N3)
N2 = Ny-n;
N3 = n;
N4 = Nx-n;

% 2. Set the source
omega = 2*pi*100; % Angular frequency of the source, unit 1/ps
T = 20; % Settling time of the source
A0 = 1; % Stable amplitude of the source
Nxs = Nx/4; % Position of first point of the source in the x axis
Nys = Ny/4; % Position of first point of the source in the y axis
% If the source is a dipole source, it will be but at (Nys,Nxs) as well
% If the source is plane wave-like
L0 = Nx/4; % Length of the source
theta = (pi/180)*30; % Oblique angle of the source, unit: radian
% Select the type of source: 1 = dipole, 2 = plane wave, 3 = Gaussian
source = 2;

% 3. % Start simulation
while(1) % Keep the simulation running until the user cancels it
    i = i+1;
    % Insert the source, which is sinusoidal but initially increases from 0
    % and gradually reaches a stable amplitude after a settling time
    if source == 1 % If the source is a dipole source 
        Hz(Nys,Nxs) = Hz(Nys,Nxs)+A0*(1-exp(-(i/T)^2))*cos(omega*i*dt);
    end
    if source == 2 % If the source is a plane wave
        for j = 0:1:L0
            Hz(Nys-round(j*sin(theta)),Nxs+round(j*cos(theta))) = ...
                Hz(Nys-round(j*sin(theta)),Nxs+round(j*cos(theta)))+...
                A0*(1-exp(-(i/T)^2))*cos(omega*i*dt);
        end
    end
    % Update the fields
    Hz(N1:N2,N3:N4) = A(N1:N2,N3:N4).*Hz(N1:N2,N3:N4)-B(N1:N2,N3:N4).*(Ey(N1:N2,N3:N4)-Ey(N1:N2,N3-1:N4-1))+...
        B(N1:N2,N3:N4).*(Ex(N1:N2,N3:N4)-Ex(N1-1:N2-1,N3:N4));
    Ex(N1:N2,N3:N4) = C(N1:N2,N3:N4).*Ex(N1:N2,N3:N4)+D(N1:N2,N3:N4).*(Hz(N1+1:N2+1,N3:N4)-Hz(N1:N2,N3:N4));
    Ey(N1:N2,N3:N4) = C(N1:N2,N3:N4).*Ey(N1:N2,N3:N4)+D(N1:N2,N3:N4).*(Hz(N1:N2,N3:N4)-Hz(N1:N2,N3+1:N4+1));
    % Display the results
    figure(1)
    image = imagesc(x,y,real(Hz),[-1 1]);
    colorbar
    colormap jet
end