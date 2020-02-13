
p1 = dlmread('dados_pressao.dat');

p = zeros(400,1024);
p(:,1) = p1(1:400,3);
for k = 2:1024
    p(:,k) = p1((k-1)*400:400*k-1,3);
end
    
t = zeros(1024,1);
t(1) = p1(1,1);
for k = 2:1024
    t(k) = p1(400*k-1,1);
end


L = 1024;
w = hann(L)';
Fs = 1/(t(2)-t(1));
T = t(2) - t(1);
X = p(200,:)-mean(p(200,:));
X = X.*w;
Yfft = fft(X);
P2 = abs(Yfft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure(1)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
[val, idx] = max(P1);

pmax = zeros(400,1);
for k = 1:400
	X = p(k,:)-mean(p(k,:));
    X = X.*w;
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    pmax(k) = Yfft(idx);
end

figure(2)
plot(pmax)

%% Dados de posição
% Posição da fonte
Y = dlmread('coordenadas_aerofolio.dat');


% Posição do observador
X = [5*cos(linspace(0,2*pi,360))',5*sin(linspace(0,2*pi,360))'];
M = 0.3;

Po = abs(intacustico(pmax,X,Y,M,f(idx)));

polarplot(linspace(0,2*pi,360),Po)


% 
% 
% n = [[X(2:end,2)-X(1:end-1,2)],[X(2:end,1)-X(1:end-1,1)]];
% n(:,1) = -n(:,1); 
% n = n./sqrt(n(:,1).^2 + n(:,2).^2);
% 
% M = 0.3;
% k = 2;
% 
% dGdx = @(xo,yo,xf,yf) (-M*k/(8*(1-M^2)^(3/2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
%     (k*(xo-xf)./((1-M^2)*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))  ).*...
%     (besselh(-1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))...
%     -besselh(1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2)));
% 
% 
% dGdy = @(xo,yo,xf,yf) (-1i/(8 *sqrt(1-M^2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
%     (k*(1-M^2)*(yf-yo)./((1-M^2)*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))).*...
%     (besselh(-1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))...
%     -besselh(1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2)));
% 
% yo = 5*sin(linspace(0,2*pi,360));
% xo = 5*cos(linspace(0,2*pi,360));
% 
% G = dGdy(xo,yo,0,0);
% 
% polarplot(linspace(0,2*pi,360),abs(G))

% t = linspace(0,1,100000)';
% t = t.^.5*pi;
% 
% c = [cos(t),sin(t)];
% 
% F = [-c(:,1),c(:,1)];
% 
% cp = [c(2:end,1) - c(1:end-1,1),c(2:end,2) - c(1:end-1,2)];
% cp = [cp;[cp(1,1)-cp(end,1),cp(1,2)-cp(end,2)]];
% 
% I = [F(:,1).*cp(:,1)+F(:,2).*cp(:,2)];
% 
% trapz(I)






