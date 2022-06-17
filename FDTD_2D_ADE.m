%% 2D FDTD simulation program for dispersive materials with PML layers
% Written by Quang Minh Dinh

clear

%% I. Set up the parameters

% 1. Basic parameters
L = 40; % Length of the simulation space, unit: um
W = 40; % Width of the simulation space, unit: um
Nx = 400; % Number of grid points in the x axis
Ny = 400; % Number of grid points in the y axis
dx = L/Nx; % Size of a pixel, unit: um
c0 = 300; % Speed of light in free-space, unit: um/ps
dt = 0.8/(c0*sqrt(2/dx^2)); % Time step that meets the CFL condition
eps0 = 8.85e-6; % Epsilon_0, unit: F/um
mu0 = 1.26; % Mu_0, unit: H/um
x = linspace(0,L,Nx); % List of coordinates in the x axis
y = linspace(0,W,Ny); % List of coordinates in the y axis
omega = 2*pi*100; % The wave frequency that we are investigating, unit: 1/ps

% 2. Create the permittivity, permeability, and conductivity profiles,
% which depend on the selected scenario
mu = ones(Ny,Nx); % Preset the permeability to be one everywhere
eps = ones(Ny,Nx); % Preset the permittivity to be one everywhere
sigma = zeros(Ny,Nx); % Preset the conductivity to be zero everywhere
scen = 4; % Chosen scenario
if scen == 1
    % Insert a near-zero-index slab of thickness Nm at ym
    Nm = Ny/4;
    ym = Ny/2;
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
    % Insert a double-negative slab of thickness Nm at ym
    Nm = Ny/4;
    ym = Ny/2;
    % Material properties
    eps_inf = 1; % Permittivity at very high frequency
    omegap = omega*1.35; % Plasma frequency, unit: 1/ps
    omega0 = omega*0.3; % Resonant frequency, unit: 1/ps
    gamma = 0; % Scattering rate, unit: 1/ps
    mu_inf = 1; % Permeability at very high frequency
    omegapm = omega*1.35; % Magnetic plasma frequency, unit: 1/ps
    omega0m = omega*0.3; % Magnetic resonant frequency, unit: 1/ps
    gammam = 0; % Magnetic scattering rate, unit: 1/ps
end
if scen == 3
    % Insert a normal dispersive dielectric slab of thickness Nm at ym
    Nm = Ny*2/5;
    ym = Ny/2;
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
if scen == 4
    % Couple radiation to surface plasmon polariton with Otto configuration
    % Insert a prism of angle thetap, permittivity epsp at Np
    thetap = (pi/180)*45; % Angle of the prism, unit: degrees
    epsp = 4; % Permittivity of the prism
    Np = Ny*2/5; % Position of the base of the prism
    for i = 0:Np-1
        eps(Np-i,round(1+i/tan(thetap)):round(Nx-i/tan(thetap))) = epsp;
    end
    % Insert a dispersive metallic layer of thickness Nm at ym
    Nm = Ny/10; % Thickness of the metallic surface
    ym = Np+10+Nm/2; % Position of the metallic slab
    eps_inf = 1; % Permittivity at very high frequency
    omegap = omega*1.73; % Plasma frequency, unit: 1/ps
    omega0 = 0; % Resonant frequency, unit: 1/ps
    gamma = 0; % Scattering rate, unit: 1/ps
    mu_inf = 1; % Permeability at very high frequency
    omegapm = 0; % Magnetic plasma frequency, unit: 1/ps
    omega0m = 0; % Magnetic resonant frequency, unit: 1/ps
    gammam = 0; % Magnetic scattering rate, unit: 1/ps
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
M = dt./(mu0*dx*(mu+1/2*dt*sigma_star));
C = eps*eps0./(eps*eps0+sigma*dt);
N = dt./(dx*(eps*eps0+sigma*dt));

% We will have the Maxwell equations like this:
% Hz = A.*Hz-M.*(Ey-Ey2)+M.*(Ex-Ex2), Ex = C.*Ex+N.*(Hz2-Hz), Ey = C.*Ey+N.*(Hz-Hz2)

% 5. Initialize the fields to to zero everywhere
Hz = zeros(Ny,Nx);
Ex = zeros(Ny,Nx);
Ey = zeros(Ny,Nx);
Bz = zeros(Ny,Nx); % Magnetic flux density field in the z axis
Jz = zeros(Ny,Nx); % Magnetic polarization density field in the z axis
Jz1 = zeros(Ny,Nx); % Jz field at one previous iteration
Jz2 = zeros(Ny,Nx); % Jz field at two previous iteration
Dx = zeros(Ny,Nx); % Displacement field in the x axis
Dy = zeros(Ny,Nx); % Displacement field in the y axis
Px = zeros(Ny,Nx); % Polarization density in the x axis
Px1 = zeros(Ny,Nx); % Px field in one previous iteration
Px2 = zeros(Ny,Nx); % Px field in two previous iteration
Py = zeros(Ny,Nx); % Polarization density field in the y axis
Py1 = zeros(Ny,Nx); % Py field in one previous iteration
Py2 = zeros(Ny,Nx); % Py field in two previous iteration

%% II. Run simulation and display the results

% 1. Basic parameters
i = 0; % Time step, count the number of iterations
n = 2; % Ignore two pixels adjacent to each boundaries
N1 = n; % The fields are only updated in the region between (N1,N3), (N1,N4), (N2,N4), (N2,N3)
N2 = Ny-n;
N3 = n;
N4 = Nx-n;
% The region of the dispersive material
NF = ym-Nm/2;
NB = ym+Nm/2;
NL = NPML+1;
NR = Nx-NPML-1;

% 2. Set the plane wave source
T = 20; % Settling time of the source
A0 = 1; % Stable amplitude of the source
L0 = Nx/4; % Length of the source
theta = (pi/180)*45; % Oblique angle of the source, unit: radian
Nxs = NPML; % Position of first point of the source in the x axis
Nys = round(NPML+L0*sin(theta)); % Position of first point of the source in the y axis

% 3. Start simulation
while(1) % Keep the simulation running until the user cancels it
    % 3.1. Update the Px, Py fields of the dispersive medium
    Px(NF:NB,NL:NR) = (2-dt^2*omega0^2)*Px1(NF:NB,NL:NR)/(1+0.5*dt*gamma)+...
        (0.5*dt*gamma-1)*Px2(NF:NB,NL:NR)/(1+0.5*dt*gamma)+...
        dt^2*eps0*omegap^2*Ex(NF:NB,NL:NR)/(1+0.5*dt*gamma);
    Py(NF:NB,NL:NR) = (2-dt^2*omega0^2)*Py1(NF:NB,NL:NR)/(1+0.5*dt*gamma)+...
        (0.5*dt*gamma-1)*Py2(NF:NB,NL:NR)/(1+0.5*dt*gamma)+...
        dt^2*eps0*omegap^2*Ey(NF:NB,NL:NR)/(1+0.5*dt*gamma);
    % 3.2. Update the Ex, Ey fields in the PMLs and the non-dispersive medium
    % In front of the dispersive medium
    Ex(N1:NF-1,N3:N4) = C(N1:NF-1,N3:N4).*Ex(N1:NF-1,N3:N4)+...
        N(N1:NF-1,N3:N4).*(Hz(N1+1:NF,N3:N4)-Hz(N1:NF-1,N3:N4));
    Ey(N1:NF-1,N3:N4) = C(N1:NF-1,N3:N4).*Ey(N1:NF-1,N3:N4)+...
        N(N1:NF-1,N3:N4).*(Hz(N1:NF-1,N3:N4)-Hz(N1:NF-1,N3+1:N4+1));
    % On the left of the medium
    Ex(NF:NB,N3:NL-1) = C(NF:NB,N3:NL-1).*Ex(NF:NB,N3:NL-1)+...
        N(NF:NB,N3:NL-1).*(Hz(NF+1:NB+1,N3:NL-1)-Hz(NF:NB,N3:NL-1));
    Ey(NF:NB,N3:NL-1) = C(NF:NB,N3:NL-1).*Ey(NF:NB,N3:NL-1)+...
        N(NF:NB,N3:NL-1).*(Hz(NF:NB,N3:NL-1)-Hz(NF:NB,N3+1:NL));
    % On the right of the medium
    Ex(NF:NB,NR+1:N4) = C(NF:NB,NR+1:N4).*Ex(NF:NB,NR+1:N4)+...
        N(NF:NB,NR+1:N4).*(Hz(NF+1:NB+1,NR+1:N4)-Hz(NF:NB,NR+1:N4));
    Ey(NF:NB,NR+1:N4) = C(NF:NB,NR+1:N4).*Ey(NF:NB,NR+1:N4)+...
        N(NF:NB,NR+1:N4).*(Hz(NF:NB,NR+1:N4)-Hz(NF:NB,NR+2:N4+1));
    % Behind the medium
    Ex(NB+1:N2,N3:N4) = C(NB+1:N2,N3:N4).*Ex(NB+1:N2,N3:N4)+...
        N(NB+1:N2,N3:N4).*(Hz(NB+2:N2+1,N3:N4)-Hz(NB+1:N2,N3:N4));
    Ey(NB+1:N2,N3:N4) = C(NB+1:N2,N3:N4).*Ey(NB+1:N2,N3:N4)+...
        N(NB+1:N2,N3:N4).*(Hz(NB+1:N2,N3:N4)-Hz(NB+1:N2,N3+1:N4+1));
    % 3.3. Update Dx, Dy, Ex, Ey in the dispersive medium
    % Update Dx, Dy, Ex, and Ey fields in the dispersive medium
    Dx(NF:NB,NL:NR) = Dx(NF:NB,NL:NR)+(dt/dx)*(Hz(NF+1:NB+1,NL:NR)-Hz(NF:NB,NL:NR));
    Dy(NF:NB,NL:NR) = Dy(NF:NB,NL:NR)+(dt/dx)*(Hz(NF:NB,NL:NR)-Hz(NF:NB,NL+1:NR+1));
    Ex(NF:NB,NL:NR) = (Dx(NF:NB,NL:NR)-Px(NF:NB,NL:NR))/(eps0*eps_inf);
    Ey(NF:NB,NL:NR) = (Dy(NF:NB,NL:NR)-Py(NF:NB,NL:NR))/(eps0*eps_inf);
    % 3.4. Update Jz field in the dispersive medium
    Jz(NF:NB,NL:NR) = (2-dt^2*omega0m^2)*Jz1(NF:NB,NL:NR)/(1+0.5*dt*gammam)+...
        (0.5*dt*gamma-1)*Jz2(NF:NB,NL:NR)/(1+0.5*dt*gammam)+...
        dt^2*mu0*omegapm^2*Hz(NF:NB,NL:NR)/(1+0.5*dt*gammam);
    % 3.5. Update Hz field in the PMLs and the non-dispersive medium
    % Infront of the dispersive medium
    Hz(N1:NF-1,N3:N4) = A(N1:NF-1,N3:N4).*Hz(N1:NF-1,N3:N4)-...
        M(N1:NF-1,N3:N4).*(Ey(N1:NF-1,N3:N4)-Ey(N1:NF-1,N3-1:N4-1))+...
        M(N1:NF-1,N3:N4).*(Ex(N1:NF-1,N3:N4)-Ex(N1-1:NF-2,N3:N4));
    % On the left of the medium
    Hz(NF:NB,N3:NL-1) = A(NF:NB,N3:NL-1).*Hz(NF:NB,N3:NL-1)-...
        M(NF:NB,N3:NL-1).*(Ey(NF:NB,N3:NL-1)-Ey(NF:NB,N3-1:NL-2))+...
        M(NF:NB,N3:NL-1).*(Ex(NF:NB,N3:NL-1)-Ex(NF-1:NB-1,N3:NL-1));
    % On the right of the medium
    Hz(NF:NB,NR+1:N4) = A(NF:NB,NR+1:N4).*Hz(NF:NB,NR+1:N4)-...
        M(NF:NB,NR+1:N4).*(Ey(NF:NB,NR+1:N4)-Ey(NF:NB,NR:N4-1))+...
        M(NF:NB,NR+1:N4).*(Ex(NF:NB,NR+1:N4)-Ex(NF-1:NB-1,NR+1:N4));
    % Behind the medium
    Hz(NB+1:N2,N3:N4) = A(NB+1:N2,N3:N4).*Hz(NB+1:N2,N3:N4)-...
        M(NB+1:N2,N3:N4).*(Ey(NB+1:N2,N3:N4)-Ey(NB+1:N2,N3-1:N4-1))+...
        M(NB+1:N2,N3:N4).*(Ex(NB+1:N2,N3:N4)-Ex(NB:N2-1,N3:N4));
    % 3.6. Insert the source, which is sinusoidal but initially increases 
    % from 0 and gradually reaches a stable amplitude after a settling time
    for j = 0:1:L0
        Hz(Nys-round(j*sin(theta)),Nxs+round(j*cos(theta))) = ...
            Hz(Nys-round(j*sin(theta)),Nxs+round(j*cos(theta)))+...
            A0*(1-exp(-(i/T)^2))*cos(omega*i*dt);
    end
    % 3.7. Update B and H fields in the dispersive medium
    Bz(NF:NB,NL:NR) = Bz(NF:NB,NL:NR)-(dt/dx)*(Ey(NF:NB,NL:NR)-Ey(NF:NB,NL-1:NR-1))+...
        (dt/dx)*(Ex(NF:NB,NL:NR)-Ex(NF-1:NB-1,NL:NR));
    Hz(NF:NB,NL:NR) = (Bz(NF:NB,NL:NR)-Jz(NF:NB,NL:NR))/(mu0*mu_inf);
    % 3.8. Save the values of Jz1, Jz2, Px1, Px2, Py1, Py2 for the next
    % iterations
    Jz2(NF:NB,NL:NR) = Jz1(NF:NB,NL:NR);
    Jz1(NF:NB,NL:NR) = Jz(NF:NB,NL:NR);
    Px2(NF:NB,NL:NR) = Px1(NF:NB,NL:NR);
    Px1(NF:NB,NL:NR) = Px(NF:NB,NL:NR);
    Py2(NF:NB,NL:NR) = Py1(NF:NB,NL:NR);
    Py1(NF:NB,NL:NR) = Py(NF:NB,NL:NR);
    % 3.9. Display results
    figure(1)
    image = imagesc(x,y,real(Hz),[-1 1]);
    colorbar
    colormap jet
    
    i = i+1;

end