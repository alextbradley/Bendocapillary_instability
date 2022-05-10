function dydx = odes(x,y, xm0, Hem, nu, k, pars)
%Returns the derviative dy/dx where 
% y = [fw,fw',fw'',...,fw'''''']= [fw,fw1,fw2,...,fw6]
% p0    :   base case pressure
% xm0   :   position of the base state meniscus
% Hem   :   channel separation of base state at meniscus
% k     :   wavenumber
% nu    :   elastocapillary number
% pars  :   parameters (= sigma)
sigma = pars(1);
    
%Equilibrium data
    He  = Hem + nu/24 /Hem *(4*xm0^3 *(xm0 - x) - (xm0 - x).^4);
    dHe = nu/24 /Hem *(-4*xm0^3 + 4*(xm0 - x).^3); %derivative on the mesh

%Extract
    fw  = y(1);
    fw1 = y(2);
    fw2 = y(3);
    fw3 = y(4);
    fw4 = y(5);
    fw5 = y(6);

%Derviatives    
    dydx = zeros(6,1);
    dydx(1) = fw1; %fw' = fw1
    dydx(2) = fw2; %fw'' = fw2
    dydx(3) = fw3;
    dydx(4) = fw4;
    dydx(5) = fw5; %etc
    dydx(6) = 1/He^3 *(3*abs(nu)*sigma*fw - 3*He^2*dHe*(fw5 - 2*k^2 *fw3 + k^4 *fw1) +...
                        -He^3*(-3*k^2*fw4 + 3*k^4 * fw2 - k^6*fw));