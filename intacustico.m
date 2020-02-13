function [Po] = intacustico(pmax,X,Y,M,omega)
%INTACUSTICO
% Resolve a integral de superfície FW-H para fonte tipo dipolo
% Entrada: 
%        p vetor das amplitudes de pressão em cada ponto numa frequencia
%        específica
%        X posição dos receptores 
%        Y posição dos pontos  na superfície.
%        M número Mach



% CONSTANTES
k = 2*pi*omega; % número de onda acústico
%% Termos da integral


% calcula normais a superficie

r = [[Y(2:end,1)-Y(1:end-1,1)],[Y(2:end,2)-Y(1:end-1,2)]];
r = [r;[Y(1,1)-Y(end,1), Y(1,2)-Y(end,2)]];
n = [[Y(2:end,2)-Y(1:end-1,2)],[Y(2:end,1)-Y(1:end-1,1)]];
n = [n;[Y(1,2)-Y(end,2),Y(1,1)-Y(end,1)]];
n(:,1) = -n(:,1); 
n = n;%./sqrt(n(:,1).^2 + n(:,2).^2);

% função de Green
%abaixo x e y são as posições da fonte

% fonte dipolar
dGdx = @(xo,yo,xf,yf) -(-M*k/(8*(1-M^2)^(3/2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
    besselh(0,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))+ ...
    (1i/(4*sqrt(1-M^2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
    (k*(xo-xf)./((1-M^2)*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))  ).* ...
    (besselh(-1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))...
    -besselh(1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2)));


dGdy = @(xo,yo,xf,yf) -(-1i/(8 *sqrt(1-M^2)))*exp(1i*M*k*(xo-xf)/(1-M^2)).* ...
    (k*(1-M^2)*(yf-yo)./((1-M^2)*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))).*...
    (besselh(-1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2))...
    -besselh(1,2,(k/(1-M^2))*sqrt((xo-xf).^2+(1-M^2)*(yo-yf).^2)));

% vetor com os resultados da função de Green 
Gx = zeros(length(Y),length(X));
Gy = zeros(length(Y),length(X));

% Gx1 = zeros(length(Y),length(X));
% Gy1 = zeros(length(Y),length(X));
% 
% Gx2 = zeros(length(Y),length(X));
% Gy2 = zeros(length(Y),length(X));
% 
% Gx3 = zeros(length(Y),length(X));
% Gy3 = zeros(length(Y),length(X));
% 

for fonte = 1:length(Y)
    for observ = 1:length(X)
        Gx(fonte,observ) = dGdx(X(observ,1),X(observ,2),Y(fonte,1),Y(fonte,2)); 
        Gy(fonte,observ) = dGdy(X(observ,1),X(observ,2),Y(fonte,1),Y(fonte,2));      
    end
end    


% x = [-sqrt(3)/5,0,sqrt(3)/5];  %pontos da quadratura
% w = [5/9,8/9,5/9]; %pesos da quadratura
% 
% aux = length(Y);
% for fonte = 1:aux
%     for observ = 1:length(X)
%         if fonte < aux
%             a = [Y(fonte,1),Y(fonte,2)];
%             b = [Y(fonte+1,1),Y(fonte+1,2)];
%         else
%             a = [Y(fonte,1),Y(fonte,2)];
%             b = [Y(1,1),Y(1,2)];            
%         end
%         
%         Gx1(fonte,observ) = dGdx(X(observ,1),X(observ,2),.5*(b(1)-a(1))*x(1) + .5*(b(1)+a(1)),.5*(b(2)-a(2))*x(1) + .5*(b(2)+a(2))); 
%         Gy1(fonte,observ) = dGdy(X(observ,1),X(observ,2),.5*(b(1)-a(1))*x(1) + .5*(b(1)+a(1)),.5*(b(2)-a(2))*x(1) + .5*(b(2)+a(2)));
%         
%         Gx2(fonte,observ) = dGdx(X(observ,1),X(observ,2),.5*(b(1)-a(1))*x(2) + .5*(b(1)+a(1)),.5*(b(2)-a(2))*x(2) + .5*(b(2)+a(2))); 
%         Gy2(fonte,observ) = dGdy(X(observ,1),X(observ,2),.5*(b(1)-a(1))*x(2) + .5*(b(1)+a(1)),.5*(b(2)-a(2))*x(2) + .5*(b(2)+a(2)));
%         
%         Gx3(fonte,observ) = dGdx(X(observ,1),X(observ,2),.5*(b(1)-a(1))*x(3) + .5*(b(1)+a(1)),.5*(b(2)-a(2))*x(3) + .5*(b(2)+a(2))); 
%         Gx3(fonte,observ) = dGdy(X(observ,1),X(observ,2),.5*(b(1)-a(1))*x(3) + .5*(b(1)+a(1)),.5*(b(2)-a(2))*x(3) + .5*(b(2)+a(2)));        
%         
%     end
%end

% dado da pressão em cada fonte multiplicado por Green
Gx = Gx.*pmax;
Gy = Gy.*pmax;

% Gx1 = Gx1.*pmax;
% Gy1 = Gy1.*pmax;
% Gx2 = Gx2.*pmax;
% Gy2 = Gy2.*pmax;
% Gx3 = Gx3.*pmax;
% Gy3 = Gy3.*pmax;


Po = zeros(1,length(X));
for observ = 1:length(X)
    %normal na direção 
%     for f = 1:length(Y)
%         if fonte < aux
%             a = [Y(f,1),Y(f,2)];
%             b = [Y(f+1,1),Y(f+1,2)];
%         else
%             a = [Y(f,1),Y(fonte,2)];
%             b = [Y(1,1),Y(1,2)];            
%         end
%         Po(observ) = Po(observ) + ...
%             .5*(b(1)-a(1))*Gx1(f,observ).*n(f,1)*w(1)+ .5*(b(2)-a(2))*Gy1(f,observ).*n(f,2)*w(1) + ...
%             .5*(b(1)-a(1))*Gx2(f,observ).*n(f,1)*w(2)+ .5*(b(2)-a(2))*Gy2(f,observ).*n(f,2)*w(2) + ...
%             .5*(b(1)-a(1))*Gx3(f,observ).*n(f,1)*w(3)+ .5*(b(2)-a(2))*Gy3(f,observ).*n(f,2)*w(3);
%     end
    I = Gx(:,observ).*n(:,1)+Gy(:,observ).*n(:,2); 
    Po(observ) = trapz(I);   
end

end

