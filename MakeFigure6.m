%Make figure 6 in the manuscript.
%(a) plot of growth rate as a function of k for nu = 5, a = 0.01 with inset
%log-log scale close to psoitive branch. 
%(b) plot of H(x0), the value of shape perturbation at the interface
%(c), (d) are as in (a), (b) but for corresponding non-wetting results,
%with nu = -5, a = -0.01.
%% Preliminaries
clear
close all
addpath('./Utilities')
figpref(4);
numplot = 1;

%parameters 
A = 1; %perturbation amplitude
V = [0.1, 0.15, 0.2,0.25, 0.3,0.35]; 
eta = 0.5;
lv = length(V);
nu = 5;
a = 0.01;
kscale = (nu*V.^3 /a).^(1/2); %scaling for when curvatures interact
sigmascale = (nu^2 * min(V).^7 ./a); %scaling for maximum growth rate
%                               (used in the accuracy of solutions)
%Search/plot parameters
krange = [logspace(-5,-2,10), 0.01:0.01:3]*max(kscale);
lk = length(krange);

%% Generate Data: Wetting case
%initialise storage
sigmaplus_positive_nu = nan(lv, lk);
h1_positive_nu = nan(lv,lk); %store value of perturbation at meniscus position

%generate data
for i = 1:length(V)
    V_ = V(i);
    
    %initialise matrix
    sigmaout = zeros(1, lk);
    
    %equilibrium details:
    [Hem, xmeq] = pleqV(V_, nu, 0); %0 suppresses plotting
    M = length(Hem);
    xms(i) = xmeq; %store the meniscus position?
    %error if more than one solution
    if M > 1
        error('found more than one equilibrium at V = %.2f, nu = %.2f', V_, nu)
    end
    
    %loop over  k values
    for q = 1:lk
        %mesh (initially) has 500 points in wet and dry
        x_mesh_wet = linspace(0,xmeq,500);
        x_mesh_dry = linspace(xmeq,1,500);
        x_mesh     = [x_mesh_wet, x_mesh_dry];
        
        %initial guess at shape and growth rate
        parsinit = 0; %initial guess at growth rate
        yinit = ones(6,1); %guess constant in each derivative
        
        %use scaling to specify accuracy
        abstol = 0.000001 * sigmascale;
        
        %solve BVP
        myodes = @(x,y,region, pars) ODEsFull(x,y,region, xmeq, Hem, nu, krange(q), pars);
        mybcs = @(yleft,yright, pars) BCsFull(yleft,yright, xmeq, Hem, nu, krange(q), a, A,eta, pars);
        options = bvpset('AbsTol', abstol, 'RelTol', 1e-5);
        solinit = bvpinit(x_mesh, yinit, parsinit); %specify mesh and guess
        solBVP = bvp4c(myodes, mybcs, solinit, options);
        
        %store data: growth rate and displacement
        sigmaout(q) = solBVP.parameters;
        shape_at_interface = deval(solBVP, xmeq-2*eps);
        h1_positive_nu(i,q) = shape_at_interface(1); %other entries are derivatives there
        
        %update progress:
        %fprintf('Completed %.f of %.f \n', count,lv*lk);
        %count = count + 1;   
    end

%add results to global info
sigmaplus_positive_nu(i,:) = sigmaout;
end

%% Non wetting case
nu = -5;
a = -0.01;
%initialise storage
sigmaplus_negative_nu = nan(lv, lk);
h1_negative_nu = nan(lv,lk); %store value of perturbation at meniscus position

%generate data
for i = 1:length(V)
    V_ = V(i);
    
    %initialise matrix
    sigmaout = zeros(1, lk);
    
    %equilibrium details:
    [Hem, xmeq] = pleqV(V_, nu, 0); %0 suppresses plotting
    M = length(Hem);
    xms(i) = xmeq; %store the meniscus position?
    %error if more than one solution
    if M > 1
        error('found more than one equilibrium at V = %.2f, nu = %.2f', V_, nu)
    end
    
    %loop over  k values
    for q = 1:lk
        %mesh (initially) has 5 points in wet and dry
        x_mesh_wet = linspace(0,xmeq,5);
        x_mesh_dry = linspace(xmeq,1,5);
        x_mesh     = [x_mesh_wet, x_mesh_dry];
        
        %initial guess at shape and growth rate
        parsinit = 0; %initial guess at growth rate
        yinit = ones(6,1); %guess constant in each derivative
        
        %use scaling to specify accuracy
        abstol = 0.000001 * sigmascale;
        
        %solve BVP
        myodes = @(x,y,region, pars) ODEsFull(x,y,region, xmeq, Hem, nu, krange(q), pars);
        mybcs = @(yleft,yright, pars) BCsFull(yleft,yright, xmeq, Hem, nu, krange(q), a, A,eta, pars);
        options = bvpset('AbsTol', abstol, 'RelTol', 1e-5);
        solinit = bvpinit(x_mesh, yinit, parsinit); %specify mesh and guess
        solBVP = bvp4c(myodes, mybcs, solinit, options);
        
        %store data: growth rate and displacement
        sigmaout(q) = solBVP.parameters;
        shape_at_interface = deval(solBVP, xmeq-2*eps);
        h1_negative_nu(i,q) = shape_at_interface(1); %other entries are derivatives there
        
        %update progress:
        %fprintf('Completed %.f of %.f \n', count,lv*lk);
        %count = count + 1;   
    end

%add results to global info
sigmaplus_negative_nu(i,:) = sigmaout;
end


%% make figure
cmap = parula(7);
%Plot the growth rate (positive branch)
figure(numplot); clf; hold on
for i = 1:lv
    plot(krange, sigmaplus_positive_nu(i,:), '-','color', cmap(i,:))
    plot(krange, sigmaplus_negative_nu(i,:), '--',...
        'color', cmap(i,:), 'HandleVisibility', 'off')

end
%tidy
ylim([-5,16]*1e-3)
%xticks(0:0.2:1)
%yticks((-4:2:4)*1e-4)
%title('growth rate of perturbations: negative branch');
xlabel('$k$', 'interpreter', 'latex', 'FontSize' ,20)
ylabel('$\sigma$', 'interpreter', 'latex', 'FontSize' ,20)
fig = gcf;
fig.Position(3:4) =  [560 640];
box on
legend({'$V = 0.1$', '$V = 0.15$', '$V = 0.2$', '$V = 0.25$', '$V = 0.3$', '$V = 0.35$'},...
    'interpreter', 'latex', 'location', 'northwest')
%add the box for inset loc
rectangle('Position', [0,-2*1e-4,0.5,2*1e-3], 'linestyle', '--', 'linewidth', 2)
numplot = numplot + 1;


% plot for logscale plot
figure(numplot); clf; hold on
for i = 1:lv
    plot(krange, sigmaplus_positive_nu(i,:), '-','color', cmap(i,:))
    plot(krange, sigmaplus_negative_nu(i,:), '--',...
        'color', cmap(i,:), 'HandleVisibility', 'off')

end
%tidy
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([1e-3, 0.5])
xlabel('$k$', 'interpreter', 'latex', 'FontSize' ,20)
ylabel('$\sigma$', 'interpreter', 'latex', 'FontSize' ,20)
fig = gcf;
fig.Position(3:4)  =  [560 320];
box on
numplot = numplot + 1;
grid on

% gap widths
figure(numplot); clf; hold on
for i = 1:lv
    plot(krange, h1_positive_nu(i,:), '-','color', cmap(i,:))
    pl = plot(krange, -h1_positive_nu(i,:), ...
        'color', (lv - i/2.5)/lv * [1,1,1], 'HandleVisibility', 'off'); %reflection in greyscale
    plot(krange, h1_negative_nu(i,:), '--',...
        'color', cmap(i,:), 'HandleVisibility', 'off')
pl.Color
end
%tidy
xlim([1e-3, 3])
xlabel('$k$', 'interpreter', 'latex', 'FontSize' ,20)
ylabel('$H(x_0)$', 'interpreter', 'latex', 'FontSize' ,20)

fig = gcf;
fig.Position(3:4)  =  [560 320];
box on
numplot = numplot + 1;
