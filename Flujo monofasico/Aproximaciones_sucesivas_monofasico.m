%% MÉTODO DE FLUJO DE POTENCIA POR APROXIMACIONES SUCESIVAS
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Datos del sistema de prueba
Sbase = 1000; Vbase = 11.4;
Zbase = ((1000*Vbase)^2)/(1000*Sbase);
Lineas = [1 2 0.1033 0.3127;2 3 0.2467 0.6051; 2 4 0.5469 0.8050];
Lineas(:,3:4)= Lineas(:,3:4)/Zbase;
Nodos = [1 0 1 0     0;    2 1 1 3000  950;
         3 1 1 2500  1970; 4 1 1 3510  2000];
Nodos(:,4:5)= Nodos(:,4:5)/Sbase;
%% Cálculos iniciales
NN = size(Nodos,1); NL = size(Lineas,1); Y = zeros(NN,NN); 
for i = 1:NL
    Ni = Lineas(i,1);  Nj = Lineas(i,2);
    ZL = Lineas(i,3) + 1j*Lineas(i,4);
    Y(Ni,Nj) = -1/ZL; Y(Nj,Ni) = -1/ZL; 
    Y(Ni,Ni) = Y(Ni,Ni) + 1/ZL;
    Y(Nj,Nj) = Y(Nj,Nj) + 1/ZL;
end
Ydg = Y(2:end,1); Ydd = Y(2:end,2:end);
IYdd = inv(Ydd);
%% Método de flujo de potencia por aproximaciones sucesivas
tmax = 100; epsilon = 1e-10; Vg = Nodos(1,3); 
Vdo = Nodos(2:end,3); Sd = Nodos(2:end,4) + 1j*Nodos(2:end,5);
for t = 1:tmax
    Vd = -IYdd*(inv(diag(conj(Vdo)))*conj(Sd) + Ydg*Vg);
    if max(abs(abs(Vd) - abs(Vdo))) < epsilon
        Vr = [Vg;Vd];
        break    
    else
    Vdo = Vd;
    end
end
%% Impresión de resultados
Ploss = real(Vr.'*conj(Y*Vr))*Sbase;
fprintf('Numero de iteraciones = %i\n',t);
fprintf('\nVoltajes en los nodos(pu)\n');
for n = 1: NN
    fprintf ('\nV%d = %.4f < %.4f \n',n,abs(Vr(n,1)),...
    angle(Vr(n,1))*180/pi)
end
fprintf('\nPerdidas de potencia en el sistema(kW)\n');
fprintf('\nPloss = %.4f \n',Ploss);