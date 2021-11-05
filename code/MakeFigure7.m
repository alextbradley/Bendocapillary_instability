%Plot the dispersion relation for different values of V.

%% Preliminaries
clear
addpath('./Utilities')
figpref(4);
% Parameters
a = 0.01;
nu = [0.01, 0.05,0.1,0.5,1,5, 10, 50, 100];
ln = length(nu);
V = [0.1, 0.2, 0.3];
lv = length(V);
lk = 100; %number of k values
A = 1; 
eta = 0.5;
colmap = parula(ln);

%storage
sigma_out = nan(lv,ln,lk);
sigmascale_out = nan(lv,ln);
krange_out = nan(lv, ln, lk);
kscale_out = nan(lv, ln);


%% Generate data
for i = 1:lv
    for j = 1:ln
        sprintf('Looking at nu = %.2f \n',nu(j))
        %problem scales
        kscale = sqrt(nu(j)*V(i)^3/a);
        kscale_out(i,j) = kscale;
        sigmascale = (nu(j)^2 * V(i).^7 ./a); %scaling for maximum growth rate
        sigmascale_out(i,j) = sigmascale;
        krange = linspace(0,kscale,lk);
        krange_out(i,j,:) = krange;
        
        %equilibrium details
        [Hem, xmeq] = pleqV(V(i), nu(j), 0); %0 suppresses plotting
      
        %setup bvp solver
        for p = 1:lk
        x_mesh_wet = linspace(0,xmeq,5);
        x_mesh_dry = linspace(xmeq,1,5);
        x_mesh     = [x_mesh_wet, x_mesh_dry];
        parsinit = 0; %initial guess at growth rate
        yinit = ones(6,1); %guess constant in each derivative
        abstol = 1e-5 * sigmascale; %use scaling to specify accuracy
        
        %solve BVP
        myodes = @(x,y,region, pars) ODEsFull(x,y,region, xmeq, Hem, nu(j), krange(p), pars);
        mybcs = @(yleft,yright, pars) BCsFull(yleft,yright, xmeq, Hem, nu(j), krange(p), a, A,eta, pars);
        options = bvpset('AbsTol', abstol, 'RelTol', 1e-5);
        solinit = bvpinit(x_mesh, yinit, parsinit); %specify mesh and guess
        solBVP = bvp4c(myodes, mybcs, solinit, options);
        
        %store sigma
        sigma_out(i,j,p) = solBVP.parameters;
        end
    end
end
      
%% Make plots
numplot = 1;
figure(numplot); clf;
for i = 1:lv
    subplot(1,lv,i); hold on
    plot([0,1], [0,0], 'k--')
    for j = 1:ln
        plot(squeeze(krange_out(i,j,:))/kscale_out(i,j), squeeze(sigma_out(i,j,:))/sigmascale_out(i,j),...
            'color', colmap(j,:));
    end

box on
xlabel('$k / k_c$', 'interpreter', 'latex', 'FontSize', 20);
ylabel('$\sigma / \sigma_c$', 'interpreter', 'latex',  'FontSize', 20);   
ylim([-0.15, 0.050001])

%add colorbar to final fig
if i == 3
    c = colorbar;
    c.Label.String = '$\nu$';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 20;
    c.Ticks = 0:0.25:1; %five ticks
    c.TickLabels = {"10^{-2}","10^{-1}","10^{0}","10^{1}","10^{2}"};
end
ax = gca;
ax.Position(2:4) = [0.16,0.22,0.75];
end
fig = gcf; fig.Position(3:4) = [1200,350];


    
        