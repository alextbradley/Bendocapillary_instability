% 30/05/22
% Convergence plot (fig A1) in Bradley et al. (2022)
%
% Plot the 2 norm of the difference between

%% Preliminaries
clear
close all
addpath('./Utilities')
numplot = 1;

set(0,'DefaultTextInterpreter','latex','DefaultAxesFontSize',14,'DefaultTextFontSize',14);
set(0, 'defaultaxesfontsize', 12);
set(0, 'defaultaxeslinewidth', 1.25);
set(0, 'defaultlinelinewidth', 1.5);
set(0, 'defaultpatchlinewidth', 0.7); 

figure(1); clf; hold on
%% Parameters
V   = [0.15,0.25,0.35]; %cross sectional volume
eta = 0.5; %Poisson's ratio
nu  = 5;   %elastocapillary number
a   = 0.01;%aspect ratio
A   = 1;   %amplitude of the perturbation

kscale     = (nu*V.^3 /a).^(1/2); %scaling for when curvatures interact
sigmascale = (nu^2 * V.^7 ./a); %scaling for maximum growth rate
N = 10*2.^(1:8); %number of grid points

%% Generate data
colmap = [1,0,0; 0,1,0;0,0,1]; 
for iV = 1:length(V)
%Get the equilibria
[Hem, xmeq] = pleqV(V(iV), nu, 0); %0 suppresses plotting
M = length(Hem);
%error if more than one solution
if M > 1
    error('found more than one equilibrium at V = %.2f, nu = %.2f', V, nu)
end

% loop over N values.
solBVP = cell(1,length(N));
errs   = zeros(1,length(N)-1);
for iN = 1:length(N)
dx = xmeq/N(iN);
N_wet = xmeq*N(iN);
N_dry = N(iN) - N_wet;

%mesh (initially) has 5 points in wet and dry
x_mesh_wet = linspace(0,xmeq, N_wet);
x_mesh_dry = linspace(xmeq, 1, N_dry);
x_mesh     = [x_mesh_wet, x_mesh_dry];


%initial guess at shape and growth rate
parsinit = 0; %initial guess at growth rate
yinit = ones(6,1); %guess constant in each derivative

%use scaling to specify accuracy
abstol = 0.000001 * sigmascale(iV);

%solve BVP
myodes = @(x,y,region, pars) ODEsFull(x,y,region, xmeq, Hem, nu, kscale(iV), pars);
mybcs = @(yleft,yright, pars) BCsFull(yleft,yright, xmeq, Hem, nu, kscale(iV), a, A,eta, pars);
options = bvpset('AbsTol', abstol, 'RelTol', 1e-5);
solinit = bvpinit(x_mesh, yinit, parsinit); %specify mesh and guess
solBVP{iN} = bvp4c(myodes, mybcs, solinit, options);

%compute the relative error against half grid size
if iN > 1
   sol_fine = solBVP{iN}.y(1,:); %solution for the shape on the current res
   sol_coarse = solBVP{iN-1}.y(1,:); %ditto for previous (half)resolution
   grid_coarse = solBVP{iN-1}.x;
   sol_fine_on_coarse = deval(solBVP{iN}, grid_coarse);  %fine solution evaluated on coarse grid. Returns all 6 shape components
   diff_coarse = sol_coarse - sol_fine_on_coarse(1,:); %difference between solutions on the coarse grid
   errs(iN - 1) = norm(diff_coarse, 2);
end


end


%% Plot the solution
plot(N(2:end), errs, 'o', 'markersize', 4, 'markeredgecolor' , colmap(iV,:));
legendinfo{iV} = ['$V = ' num2str(V(iV)), '$'];
end

%% Tidy plot
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on; box on
xlabel('$N$');
ylabel('$||\mathbf{Y}_{N} - \mathbf{Y}_{N/2}||_{2}$')
legend(legendinfo, 'interpreter', 'latex')

