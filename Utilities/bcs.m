function res = bcs(ya,yb, xm0, Hem, nu, k, a, A,eta, pars)
%Return the residuals from the boundary conditions. 
% ya    :   solution (y = [fw,...,fw''''']) at x = a = 0
% yb    :   solution (y = [fw,...,fw''''']) at x = b = xm0
% xm0   :   base case meniscus position
% Hem   :   base case meniscus gap
% nu    :   elastocapillary number
% k     :   wavenumber
% a     :   (signed) aspect ratio
% A     :   amplitude of perturbation
% eta   :   Poisson's ratio
% pars  :   parameters (= sigma)

sigma = pars(1);

    
res = zeros(7,1); %init
%equilibrium data
    p0 = -nu/Hem;
    dHem = -nu/24 /Hem * 4 *xm0^3; %slope at x = xm0

%BC at x = 0:
    fw   = ya(1);
    dfw  = ya(2);
    d2fw = ya(3);
    d3fw = ya(4);
    d4fw = ya(5);
    d5fw = ya(6);
    
    %clamping
    res(1) = fw;
    res(2) = dfw;   
    
    %no flux condition
    res(3) = d5fw - 2 * k^2 * d3fw + k^4 * dfw; 
%BC at x = xm0 = b
    fw = yb(1);
    dfw  = yb(2);
    d2fw = yb(3);
    d3fw = yb(4);
    d4fw = yb(5);
    d5fw = yb(6);
    
    %effective bc
    res45 = bc_matching(fw, dfw, d2fw, d3fw, k, xm0, p0, A, eta);
    res(4) = res45(1);
    res(5) = res45(2);
        
    %Pressure bc
    res(6) = (d4fw - 2 * k^2 *d2fw + k^4 * fw) - nu/Hem^2 *(A*dHem + fw) - nu*a*k^2*A;
    %kinematic BC
    res(7) = 3*abs(nu)*sigma*A + Hem^2 *(d5fw - 2*k^2*d3fw + k^4*dfw);
    %d5fw - 2*k^2*d3fw + k^4*dfw
        
    
     
      