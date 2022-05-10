function res = BCsFull(yleft,yright, xm0, Hem, nu, k, a, A,eta, pars)
%Return the residuals from the boundary conditions.
% xm0   :   base case meniscus position
% Hem   :   base case meniscus gap
% nu    :   elastocapillary number
% k     :   wavenumber
% a     :   (signed) aspect ratio
% A     :   amplitude of perturbation
% eta   :   Poisson's ratio
% pars  :   parameters (= sigma)

%yleft(:,1)  : wet solution at x = 0
%yright(:,1) : wet solution at x = xm
%yleft(:,2)  : dry solution at x = xm
%yright(:,2) : dry solution at x = 1;
sigma = pars(1);
res = zeros(11,1);

%clamping
res(1) = yleft(1,1); %f = 0
res(2) = yleft(2,1); %f' = 0

%no flux at x = 0
res(3) = yleft(6,1) - 2*k^2 * yleft(4,1) + k^4*yleft(2,1);

%continuity
res(4) = yright(1,1) - yleft(1,2); %continuity of f across x = xm
res(5) = yright(2,1) - yleft(2,2); %continuity of f' across x = xm 
res(6) = yright(3,1) - yleft(3,2); %continuity of f'' across x = xm 
res(7) = yright(4,1) - yleft(4,2) - A*nu/Hem; %jump in f''' across x = xm
% fw''' - fd''' + A*p0 = 0 where p0 = -nu/Hem is base case pressure

%pressure condition
dHem = -nu/24 /Hem * 4 *xm0^3; %slope at x = xm0
res(8) = (yright(5,1) - 2 * k^2 *yright(3,1) + k^4 * yright(1,1))...
    - nu/Hem^2 *(A*dHem + yright(1,1)) - nu*a*k^2*A;

%kinematic bc
res(9) = 3*abs(nu)*sigma*A + Hem^2 *(yright(6,1) - 2*k^2*yright(4,1) + k^4*yright(2,1));

%free ends at x = 1
res(10) = yright(3,2) - eta*k^2*yright(1,2);
res(11) = yright(4,2) - (2-eta)*k^2*yright(2,2);

%artificial conditions f''''' = 0, f'''''' = 0 at x = 1
res(12) = yright(5,2);
res(13) = yright(6,2);
