%Make figure 3c and 3d (equilibria characterization and growth rate for
%non-wetting case)

%% Parameters
nu = -5;
%guess at sigma for V = 1
sigmaguess = -43.4; %nu = -5

Vrange = 0.40:0.5*1e-4:1; %need fine resolution on V to facilitate continuation method
%make sure Vrange has 1 as upper limit(can shrink later if you like)
%need very fine resolution on V for continuation
sigmaout = nan(1,length(Vrange)); %only one solution for non-wetting

%other parameters
a = 0.01; %aspect ratio
eta = 0.4; %poisson's ratio
k = 0; %uniform perturbation only
A = 1; %amplitude of perturbation (=1 wlog)
plotit = 0;

%% get equilibrium info
numsols = nan(1,length(Vrange));
xms = nan(1,length(Vrange));
hms = nan(1,length(Vrange));
h1s = nan(1,length(Vrange));

for i = 1:length(Vrange)
    [hm, xm] = pleqV(Vrange(i), nu, plotit);
    numsols(i) = length(hm); 
    xms(i) = xm;
    hms(i) = hm;
    h1s(i) = hm + nu*xm.^3 /6 ./hm .*(xm -1); %channel width at x = 1
end %end loop over V

%% Get the solution at V = 1
[~,idx] = min(abs(Vrange - 1));
xmeq = xms(idx);
Hem = hms(idx);

%initialise

parsinit = sigmaguess; %initial guess at growth rate
xmesh = linspace(0,xmeq,50);
omega = 2*pi*5/xmeq; %initial guess is periodic with 5 periods (this just seems to work well as an initial guess)
guess = @(x) mat4init(x, omega); %initial guess is periodic

%Solve
myodes = @(x,y,pars) odes(x,y, xmeq, Hem, nu, k, pars);
mybcs = @(y0,yxm0, pars) bcs(y0,yxm0, xmeq, Hem, nu, k, a, A,eta, pars);
solinit = bvpinit(xmesh, guess, parsinit); %specify mesh and guess
sol = bvp4c(myodes, mybcs, solinit);


%% extend solution down to smaller V
for i = flip(1:length(Vrange))
   Hem = hms(i);
   xmeq = xms(i);
    
   % get initial guess 
   if i == length(Vrange) %if we're adjacent to two sols region
       solBVP = sol; 
   end 
   
   %solve the bvp
   parsinit = solBVP.parameters; %initial guess at growth rate
   xmesh = linspace(0,xmeq,50);
   omega = 2*pi*5/xmeq; %five periods
   guess = @(x) mat4init(x, omega); %initial guess is periodic
   solinit = bvpinit(xmesh, guess, parsinit); %specify mesh and guess
   myodes = @(x,y,pars) odes(x,y, xmeq, Hem, nu, k, pars);
   mybcs = @(y0,yxm0, pars) bcs(y0,yxm0, xmeq, Hem, nu, k, a, A,eta, pars);
   solBVP = bvp4c(myodes, mybcs, solinit);
   
   %store the solution
   sigmaout(i) = solBVP.parameters;
   Vrange(i)
end
%% plots
clf;
colmap = [0, 47, 167; 
         167, 0, 47]/255;

 i = 2; %use 'stable' color    
subplot(1,2,1); hold on
plot(Vrange, h1s, 'color', colmap(i,:), 'linewidth', 2.5)
box on
xlabel('$V$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$h_e(x=1)$', 'interpreter', 'latex', 'FontSize', 20)
xlim([0.45, 0.7+1e-7])

   
subplot(1,2,2); hold on
idx = sigmaout < -1; %don't plot cases where we pick up zero root
plot(Vrange(idx), sigmaout(idx), 'color', colmap(i,:), 'linewidth', 2.5)
%add zero root
plot(Vrange, zeros(1, length(Vrange)), 'color', colmap(i,:), 'linewidth', 2.5)  
box on
xlabel('$V$', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\sigma_u$', 'interpreter', 'latex', 'FontSize', 20)
xlim([0.45, 0.7+1e-7])
ylim([-800, 100])
