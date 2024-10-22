%% MÉTODO DE FLUJO DE POTENCIA POR BARRIDO ITERATIVO
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Datos del sistema de prueba
Sbase = 1000; Vbase = 13.2; 
Zbase = ((1000*Vbase)^2)/(1000*Sbase);
Lineas = [1 2 0.3968 0.5290; 1 4 0.4232 0.5819; 2 3 0.4761 0.9522;
          3 4 0.2116 0.2116; 3 5 0.4232 0.5819; 4 5 0.5819 0.8050];
Lineas(:,3:4) = Lineas(:,3:4)/Zbase;
Nodos = [1 0 1 0    0; 2 1 1 2000 1250; 3 1 1 2500 1550; 
         4 1 1 3000 1900; 5 1 1 3800 2310];
Nodos(:,4:5) = Nodos(:,4:5)/Sbase;
%% Cálculos iniciales
NN = size(Nodos,1); NL = size(Lineas,1); A = zeros(NN,NL);
Zp = zeros(NL,NL);
for i = 1:NL
    Ni = Lineas(i,1); Nj = Lineas(i,2);
    A(Ni,i)= 1; A(Nj,i)= -1;
    Zp(i,i) = Lineas(i,3) + 1j*Lineas(i,4);
end
Ag = A(1,:); Ad = A(2:end,:); 
Yp = inv(Zp); Zdd = inv(Ad*Yp*Ad.');
%% Método de flujo de potencia por barrido iterativo
tmax = 100; epsilon = 1e-10; Vg = Nodos(1,3); 
Vdo = Nodos(2:end,3); Sd = Nodos(2:end,4) + 1j*Nodos(2:end,5);
for t = 1:tmax
    Vd = -Zdd*(inv(diag(conj(Vdo)))*conj(Sd) + (Ad*Yp*Ag.')*Vg);
    if max(abs(abs(Vd) - abs(Vdo))) < epsilon
        Vr = [Vg;Vd];
        break    
    else
    Vdo = Vd;
    end
end
%% Impresión de resultados
E = A.'*Vr; J = Yp*E;
Ploss = real(J.'*Zp*conj(J))*Sbase;
fprintf('Número de iteraciones = %i\n',t);
fprintf('\nVoltajes en los nodos(pu)\n');
for n = 1: NN
    fprintf ('\nV%d = %.4f < %.4f \n',n,abs(Vr(n,1)), ...
    angle(Vr(n,1))*180/pi)
end
fprintf('\nPérdidas de potencia en el sistema(kW)\n');
fprintf('\nPloss = %.4f \n',Ploss);