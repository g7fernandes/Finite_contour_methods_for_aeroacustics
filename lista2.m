M = 0.3;
k = 2;

dGdx = @(xo,yo,xf,yf) -(-M*k/(8*(1-M^2)^(3/2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
    besselh(0,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))+ ...
    -(1i/(4*sqrt(1-M^2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
    (k*(xo-xf)./((1-M^2)*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))  ).* ...
    (besselh(-1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))...
    -besselh(1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2)));


dGdy = @(xo,yo,xf,yf) -(-1i/(8 *sqrt(1-M^2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
    (k*(1-M^2)*(yf-yo)./((1-M^2)*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))).*...
    (besselh(-1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))...
    -besselh(1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2)));

yo = 5*sin(linspace(0,2*pi,360));
xo = 5*cos(linspace(0,2*pi,360));

G = dGdy(xo,yo,0,0);

polarplot(linspace(0,2*pi,360),abs(G))