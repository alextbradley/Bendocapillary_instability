function dydx = ODEsFull(x,y, region, xm0, Hem, nu, k, pars)
%Returns the derviative dy/dx where
% y = [fw,fw',fw'',...,fw''''',fd,fd',fd'',fd''']= [fw,fw1,fw2,...,fw6, fd,fd1,fd2,fd3]]
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
f  = y(1);
f1 = y(2);
f2 = y(3);
f3 = y(4);
f4 = y(5);
f5 = y(6);

dydx = zeros(6,1);

switch region
    case 1 % x in [0, xm]
        dydx(1) = f1; %fw' = fw1
        dydx(2) = f2; %fw'' = fw2
        dydx(3) = f3;
        dydx(4) = f4;
        dydx(5) = f5; %etc
        dydx(6) = 1/He^3 *(3*abs(nu)*sigma*f - 3*He^2*dHe*(f5 - 2*k^2 *f3 + k^4 *f1) +...
            -He^3*(-3*k^2*f4 + 3*k^4 * f2 - k^6*f));
    case 2 % x in [xm,1]
        dydx(1) = f1;
        dydx(2) = f2;
        dydx(3) = f3;
        dydx(4) = 2*k^2*f2 - k^4*f;
        %other two components are zero
end
end