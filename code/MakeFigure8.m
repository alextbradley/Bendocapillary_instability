%Make figure 8 in the manuscript:
%(a) Normalised channel width perturbation at x_0 as a function of nu*V^5/a
%(b) Normalised growth rate as a function of normalized wavenumber.
%Coloured by nu*V^5 /a, with inset of normalized sigma^* versus nu*V^5/a

%% Preliminaries
clear
addpath('./Utilities')
figpref(4);
% Parameters
a = 0.01;
nu = [0.01, 0.05,0.1,0.5,1,5, 10, 50, 100];
ln = length(nu);
V = [0.1, 0.2, 0.3, 0.4];
lv = length(V);
lk = 1e2; %number of k values (nb will take O(mins) for lk = 1e3)
A = 1; 
eta = 0.5;

%storage
sigma_out = nan(lv*ln,lk);
sigmascale_out = nan(lv*ln,1);
krange_out = nan(lv*ln, lk);
H_x0_out = nan(lv*ln, lk);
kscale_out = nan(lv*ln,1);
delta_out = nan(lv*ln,1);
eps_out = nan(lv*ln,1);
nuV3 = nan(lv*ln,1);

%colormap
nc = 1e3;
maxval = (max(nu)*max(V)^5/a);
minval = (min(nu)*min(V)^5/a);
logmaxval = log(maxval);
logminval = log(minval);
%colour by eps
maxval = (max(nu)*max(V)^4);
minval = (min(nu)*min(V)^4);
logmaxval = log(maxval);
logminval = log(minval);

colmap = parula(nc);
%% Generate data
count = 1;
for i = 1:lv
    for j = 1:ln
        nu(j)
        %problem scales
        kscale = sqrt(nu(j)*V(i)^3/a);
        kscale_out(count) = kscale;
        sigmascale = (nu(j)^2 * V(i)^7 /a); %scaling for maximum growth rate
        sigmascale_out(count) = sigmascale;
        krange = linspace(0,kscale,lk);
        krange_out(count,:) = krange;
        delta_out(count) = nu(j)*V(i)^5/a;
        eps_out(count) = nu(j)*V(i)^4;
        nuV3(count) = nu(j)*V(i)^3;
        
        %equilibrium details
        [Hem, xmeq] = pleqV(V(i), nu(j), 0); %0 suppresses plotting
        if ~isempty(Hem)
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
        
        %store sigma and channel width at x = 1
        sigma_out(count,p) = solBVP.parameters;
        [~,idx] = min(abs(solBVP.x - xmeq));
        He = solBVP.y;
        H_x0_out(count,p) = He(1,idx);
        end
        end
        count = count + 1;
    end
end
      
%% Make plots
numplot = 1;
figure(numplot); clf; 
subplot(1,2,1); hold on
for i = 1:ln*lv
    %work out the color
    %ratio = (log(delta_out(i)) - logminval)/(logmaxval - logminval); %how close to min or max value are we
    ratio = (log(eps_out(i)) - logminval)/(logmaxval - logminval); %how close to min or max value are we
    [~,idx] = min(abs(linspace(0,1,nc) - ratio));
    plot(krange_out(i,:)/kscale_out(i), 3*H_x0_out(i,:)/nuV3(i), 'color', colmap(idx, :))
end
plot([0,1], -[1,1], 'k--')
box on
xlabel('$k / k_c$', 'interpreter', 'latex', 'FontSize', 20);
ylabel('$6 H(x_0)/(\nu V^3)$', 'interpreter', 'latex',  'FontSize', 20); 
ax2 = gca; 
ax2.Position(3) = 0.25;
yticks([-2:0.5:0])


subplot(1,2,2); hold on
for i = 1:ln*lv
    %work out the color
    ratio = (log(delta_out(i)) - logminval)/(logmaxval - logminval); %how close to min or max value are we
    ratio = (log(eps_out(i)) - logminval)/(logmaxval - logminval); %how close to min or max value are we
    [~,idx] = min(abs(linspace(0,1,nc) - ratio));
    plot(krange_out(i,:)/kscale_out(i), sigma_out(i,:)/sigmascale_out(i), 'color', colmap(idx, :))
end
box on
xx = linspace(0,1);
yy = -(xx.^2.*(2*xx.^2 - 1))/6;
plot(xx, yy, 'k--')
xlabel('$k / k_c$', 'interpreter', 'latex', 'FontSize', 20);
ylabel('$\sigma/\sigma_c$', 'interpreter', 'latex',  'FontSize', 20); 
c = colorbar;
c.Label.String = '$\nu V^4$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 20;
%c.Ticks = 1/7:2/7:1; %five ticks
%c.TickLabels = {"10^{-4}","10^{-2}","10^{0}","10^{2}","10^{2}"};
c.Ticks = linspace(0,1,8);
c.TickLabels = {"10^{-5}","10^{-4}","10^{-3}","10^{-2}","10^{-1}","10^{0}","10^{1}","10^{2}"};
ax2=gca; 
ax2.Position(1) = 0.46;
ax2.Position(3) = 0.44;

axnew = axes;
axnew.Position = [0.51,0.22, 0.3, 0.33];
colmap2 = [0, 47, 167;  167, 47, 167;47, 167, 0; 167,0, 47]/255;
hold on
for i = 1:lv
    VV = V(i);
    %get the largest growth rate for each different nu at this V
    sigmamax = nan(ln,1);
    for j = 1:ln
        sigmas = sigma_out((i-1)*ln + j,:);
        sigmamax(j) = 48*max(sigmas)/sigmascale_out((i-1)*ln + j);
    end
    plot(eps_out((i-1)*ln + 1:(i-1)*ln + ln), sigmamax, 'o-','color', colmap2(i,:), 'markersize', 5)
end
xlim([1e-5, 10])
ylim([0.2, 1.1])
ax = gca;
plot(ax.XLim, [1,1], 'k--');
set(gca, 'XScale', 'log')
box on
xlabel('$\nu V^4 $', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$48 \sigma^*/\sigma^*_{SD}$', 'interpreter', 'latex',  'FontSize', 16); 

fig = gcf;
fig.Position(3:4) = [1100, 500];