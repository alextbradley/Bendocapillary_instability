function hm = pleq(xm0, nu, plotit)
%Return the equilibrium configuration data and plot the solution.

%solve the quadratic for Hm
C = [1, -1, nu*xm0^4/8];
hm = roots(C);

%remove anything unwanted
hm = hm(hm>0);
hm = hm(hm<1);
hm = hm(hm == real(hm));

%check regime 1
idx = hm + nu*xm0^3/6 ./hm *(xm0-1) >0;
hm = hm(idx);

%plot each solution (may be multiple)
if ~isempty(hm)
    for i = 1:length(hm)

        %anonymous function of wet shape
        hw = @(x)( hm(i) + nu/24 /hm(i) *(4*xm0^3*(xm0 - x) - (xm0 - x).^4));
        
        %anonymous function of dry shape
        hd = @(x) (hm(i) + nu*xm0^3 /6 /hm(i) *(xm0 -x));
        
        xw = linspace(0,xm0);
        xd = linspace(xm0,1);
        if plotit
            figure(i); clf; hold on
            plot(hw(xw),xw, 'k', hd(xd),xd, 'k')
            plot(-hw(xw),xw, 'k', -hd(xd),xd, 'k')

            %add drop
            fill([hw(xw) -flip(hw(xw))], [xw, flip(xw)], 'b')
        end
        %associated volume
        V = xm0*hm(i) + 3*nu*xm0^5 /40 /hm(i);
    end
else
    warning('No equilibria found')
end
