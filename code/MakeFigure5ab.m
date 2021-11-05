%Make figure 5a and 5b (equilibria characterization and growth rate for
%wetting case)
addpath('Utilities')
figpref(4);
%% Parameters
nu = 5;
%sigmaguess = [20.6, -428]; %guess at sigma (unstable first), at point
%where two roots first appear
%above for nu = 2
sigmaguess = [18.4, -177.4]; %nu = 5;
%sigmaguess = [18.6, -94]; %nu = 10;

Vrange = 0.45:1e-3:0.7; %need fine resolution on V to facilitate continuation method
%above values for initial guesses based on 1e-3 resolution!
sigmaout = nan(2,length(Vrange)); 

%other parameters
a = 0.01; %aspect ratio
eta = 0.4; %poisson's ratio
k = 0; %uniform perturbation only
A = 1; %amplitude of perturbation (=1 wlog)
plotit = 0;

%% find number of solutions for each value of V
numsols = nan(1,length(Vrange));
xms = nan(2,length(Vrange));
hms = nan(2,length(Vrange));
h1s = nan(2,length(Vrange));
for i = 1:length(Vrange)
    [hm, xm] = pleqV(Vrange(i), nu, plotit);
    h1 = hm + nu*xm.^3 /6 ./hm .*(xm -1); %channel width at x = 1
    numsols(i) = length(hm); 
    if length(hm) == 1
        hms(2,i) = hm;
        xms(2,i) = xm;
        h1s(2,i) = h1;
    elseif length(hm) == 2
        [hm, I] = sort(hm);
        xm = xm(I);
        hms(:,i) = hm;
        xms(:,i) = xm;
        h1s(:,i) = h1;
    end
    
end %end loop over V

%% start by getting bvp solutions at the two solution point
idx1 = find(numsols == 2, 1, 'first');
idx2 = find(numsols == 2, 1, 'last');

for i=1:2
xm = xms(i, idx1);
hm = hms(i,idx1);
%initialise
parsinit = sigmaguess(i); %initial guess at growth rate
xmesh = linspace(0,xm,50);
omega = 2*pi*5/xm; %initial guess is periodic with 5 periods (this just seems to work well as an initial guess)
guess = @(x) mat4init(x, omega); %initial guess is periodic

%Solve
myodes = @(x,y,pars) odes(x,y, xm, hm, nu, k, pars);
mybcs = @(y0,yxm0, pars) bcs(y0,yxm0, xm, hm, nu, k, a, A,eta, pars);
solinit = bvpinit(xmesh, guess, parsinit); %specify mesh and guess
sol = bvp4c(myodes, mybcs, solinit);
twoEq_bvpSol{i} = sol;
twoEq_growthrate(i) = sol.parameters;
end

%% Extend to smaller V
disp('performing one root growth rate search');
V_oneEq = Vrange(1:idx1-1);

for i=flip(1:length(V_oneEq))
   %load equilibrium details
   hm = hms(2,i);
   xm = xms(2,i);
    
   % get initial guess 
   if i == length(V_oneEq) %if we're adjacent to two sols region
       solBVP = twoEq_bvpSol{2}; 
   end
   
   %solve the bvp
   parsinit = solBVP.parameters; %initial guess at growth rate
   xmesh = linspace(0,xm,50);
   omega = 2*pi*5/xm; %five periods
   guess = @(x) mat4init(x, omega); %initial guess is periodic
   solinit = bvpinit(xmesh, guess, parsinit); %specify mesh and guess
   myodes = @(x,y,pars) odes(x,y, xm, hm, nu, k, pars);
   mybcs = @(y0,yxm0, pars) bcs(y0,yxm0, xm, hm, nu, k, a, A,eta, pars);
   solBVP = bvp4c(myodes, mybcs, solinit);
   
   %store the solution
   sigmaout(2,i) = solBVP.parameters;   
   
end


%% Extend to larger V
V_twoEq = Vrange(idx1:idx2);
disp('one root growth rate search finished, now looking at two equilibrium part');
for m = 1:2
for i = idx1:idx2
   %load equilibrium details
   hm = hms(m,i);
   xm = xms(m,i);
    
   % get initial guess 
   if i == idx1 %if we're adjacent to two sols region
       solBVP = twoEq_bvpSol{m}; 
   end
   
   %solve the bvp
   parsinit = solBVP.parameters; %initial guess at growth rate
   xmesh = linspace(0,xm,50);
   omega = 2*pi*5/xm; %five periods
   guess = @(x) mat4init(x, omega); %initial guess is periodic
   solinit = bvpinit(xmesh, guess, parsinit); %specify mesh and guess
   myodes = @(x,y,pars) odes(x,y, xm, hm, nu, k, pars);
   mybcs = @(y0,yxm0, pars) bcs(y0,yxm0, xm, hm, nu, k, a, A,eta, pars);
   solBVP = bvp4c(myodes, mybcs, solinit);
   
   %store the solution
   sigmaout(m,i) = solBVP.parameters;
   solBVP.parameters;
   
end
end

%% plots
clf;
colmap = [0, 47, 167; 
         167, 0, 47]/255;

     
subplot(1,2,1); hold on
for i = 1:2
plot(Vrange, h1s(i,:), 'color', colmap(i,:), 'linewidth', 2.5)
end
box on
xlabel('$V$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$h_e(x=1)$', 'interpreter', 'latex', 'FontSize', 20)
xlim([0.45, 0.7+1e-7])
     
subplot(1,2,2); hold on
for i = 1:2
plot(Vrange, sigmaout(i,:), 'color', colmap(i,:), 'linewidth', 2.5)
%add zero root
if i == 1
plot(Vrange(idx1:idx2), 2.5*(-1)^(i+1) * ones(1,idx2 - idx1 + 1), 'color', colmap(i,:), 'linewidth', 2)
else
 plot(Vrange(1:idx2), 2.5*(-1)^(i+1) * ones(1,idx2), 'color', colmap(i,:), 'linewidth', 2)
end   
end
%box in main fig for zoom
rectangle('Position', [0.53, -45,0.17,80], 'linestyle', '--', 'linewidth', 2)

%tidy
box on
xlim([min(Vrange), max(Vrange)]);
xlabel('$V$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\sigma_u$', 'interpreter', 'latex', 'FontSize', 20)
ylim([-800, 100])
xlim([0.45, 0.70000001])


%inset
axes('Position', [.67,.23,.22,.31]);
box on
hold on
for i = 1:2
plot(Vrange, sigmaout(i,:), 'color', colmap(i,:), 'linewidth', 2.5)
%add zero root
if i == 1
plot(Vrange(idx1:idx2), 0.5*(-1)^(i+1) * ones(1,idx2 - idx1 + 1), 'color', colmap(i,:), 'linewidth', 2)
else
 plot(Vrange(1:idx2), 0.5*(-1)^(i+1) * ones(1,idx2), 'color', colmap(i,:), 'linewidth', 2)
end
end
ylim([-40, 25])
xlim([0.53, 0.68])
shg

fig = gcf;
fig.Position(3:4) = [785 308];


