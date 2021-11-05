function [hm, xm] = pleqV(V, nu, plotit)
%returns the meniscus displacement and meniscus position associated with
%regime 1 equilibrium with volume V, bendability nu. Plotit specifies
%whether to plot or not.

%Compute roots
%Eq meniscus positions satisfy nu*xm^6 + 30*xm^2 -80*V*xm + 50*V^2 = 0
c = [nu, 0,0,0,30,-80*V, 50*V^2];
xm = roots(c); %roots of the polynomial

%remove any unwanted roots
xm = xm(real(xm)==xm);
xm = xm(xm>0);
xm = xm(xm<1);

%compute associated meniscus displacements: 
hm = (5*V - 3*xm)./(2*xm);


%Plot
if ~isempty(hm)
    %Remove any roots which violate h(x=1) > 0
    idx = ((hm + nu*xm.^3.*(xm-1)./hm./6) >0);
    hm = hm(idx);
    xm = xm(idx);
    for i = 1:length(hm)
        %anonymous function of wet shape
        hw = @(x)( hm(i) + nu/24 /hm(i) *(4*xm(i)^3*(xm(i) - x) - (xm(i) - x).^4));
        %anonymous function of dry shape
        hd = @(x) (hm(i) + nu*xm(i)^3 /6 /hm(i) *(xm(i) -x));
        %independent variables in wet and dry
        xw = linspace(0,xm(i));
        xd = linspace(xm(i),1);
        if plotit
            figure(i); clf; hold on
            plot(hw(xw),xw, 'k', hd(xd),xd, 'k')
            plot(-hw(xw),xw, 'k', -hd(xd),xd, 'k')

            %add drop
            fill([hw(xw) -flip(hw(xw))], [xw, flip(xw)], 'b')
            title(['Equilibrium ' num2str(i)]);
        end
        %associated volume
        %V = xm(i)*hm(i) + 3*nu*xm(i)^5 /40 /hm(i)
    end
else
    %warning('No equilibria found')
end