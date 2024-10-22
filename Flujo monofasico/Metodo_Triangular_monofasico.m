%% MÉTODO DE FLUJO DE POTENCIA POR MÉTODO TRIANGULAR
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Datos del sistema de prueba
Sbase = 1000; Vbase = 11.4; Zbase = ((1000*Vbase)^2)/(1000*Sbase);  
Lineas = [1 2 0.5025 0.3025; 2 3 0.4020 0.2510;
          2 4 0.3660 0.1864; 2 5 0.2872 0.4088];
Lineas(:,3:4) = Lineas(:,3:4)/Zbase;
Nodos = [1 0 1 0    0; 2 1 1 900  500; 3 1 1 1200 950; 
         4 1 1 2000 1150; 5 1 1 2500 1200];
Nodos(:,4:5) = Nodos(:,4:5)/Sbase;
%% Cálculos iniciales
NN = size(Nodos,1); NL = size(Lineas,1);
T = zeros(NL,NN); Zp = zeros(NL,NL);
for i = NL:-1:1 
    Ni = Lineas(i,1); Nj = Lineas(i,2); T(i,Nj) = 1;
    for j = NL-1:-1:1
        if Lineas(j,2) == Ni
           T(j,Nj) = 1;
           Ni = Lineas(j,1); 
        end
    end
    Zp(i,i) = Lineas(i,3) + 1j*Lineas(i,4);
end
T(:,1) = [];
Zbb = T.'*(Zp*T);
%% Método de flujo de potencia por metodo triangular
tmax = 100; epsilon = 1e-10; Vg = Nodos(1,3); 
Vdo = Nodos(2:end,3); unos = ones(size(Vdo,1),1); 
Sd = Nodos(2:end,4) + 1j*Nodos(2:end,5);
for t = 1:tmax
    Vd = unos*Vg - Zbb*(inv(diag(conj(Vdo)))*conj(Sd));
    if max(abs(abs(Vd) - abs(Vdo))) < epsilon
        Vr = [Vg;Vd];
        break    
    else
    Vdo = Vd;
    end
end
%% Impresión de resultados
Id = inv(diag(conj(Vd)))*conj(Sd); J = T*Id;
Ploss = real(J.'*Zp*conj(J))*Sbase;
fprintf('Número de iteraciones = %i\n',t);
fprintf('\nVoltajes en los nodos(pu)\n');
for n = 1: NN
    fprintf ('\nV%d = %.4f < %.4f \n',n,abs(Vr(n,1)),...
    angle(Vr(n,1))*180/pi)
end
fprintf('\nPérdidas de potencia en el sistema(kW)\n');
fprintf('\nPloss = %.4f \n',Ploss);