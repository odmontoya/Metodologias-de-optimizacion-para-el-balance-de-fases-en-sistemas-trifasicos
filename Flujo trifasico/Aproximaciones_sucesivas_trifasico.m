%% MÉTODO DE FLUJO DE POTENCIA POR APROXIMACIONES SUCESIVAS
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Datos del sistema de prueba
clc, clear
Sbase = 1000; Vbase = 7.2; Zbase = ((1000*Vbase)^2)/(1000*Sbase);
Lineas = [1 2 1; 1 3 1; 2 3 1];
C1 = exp(1j*-2*pi/3); C2 = exp(1j*2*pi/3);
Nodos = [1 0 1 C1 C2 0 0 0 0 0 0 0;
         2 1 1 C1 C2 975 605 1525 945 1300 805 0;
         3 1 1 C1 C2 1800 1100 300 185 850 526 0];           
Nodos(:,6:11) = Nodos(:,6:11)/Sbase;       
%% Cálculos iniciales
NN = size(Nodos,1); NL = size(Lineas,1); Y3 = zeros(3*NN,3*NN); 
for i = 1:NL
    Ni = Lineas(i,1); Nj = Lineas(i,2);
    ZL = MatrizZ22(Lineas(i,3))/Zbase;
    Y3(3*Ni-2:3*Ni,3*Nj-2:3*Nj) = -inv(ZL); 
    Y3(3*Nj-2:3*Nj,3*Ni-2:3*Ni) = -inv(ZL); 
    Y3(3*Ni-2:3*Ni,3*Ni-2:3*Ni)=Y3(3*Ni-2:3*Ni,3*Ni-2:3*Ni)+inv(ZL);
    Y3(3*Nj-2:3*Nj,3*Nj-2:3*Nj)=Y3(3*Nj-2:3*Nj,3*Nj-2:3*Nj)+inv(ZL);
end
Ydg3 = Y3(4:end,1:3); Ydd3 = Y3(4:end,4:end);
IYdd3 = inv(Ydd3);
%% Método de flujo de potencia por aproximaciones sucesivas
tmax = 100; epsilon = 1e-10;
Vg3 = [Nodos(1,3),Nodos(1,4),Nodos(1,5)].'; 
Vd3o = zeros(3*NN-3,1); Sd3 = zeros(3*NN-3,1);
Id3 = zeros(3*NN-3,1);
H = [0 0 1; 1 0 0; 0 1 0]; M = [1 -1 0; 0 1 -1; -1 0 1];
for k = 1:NN-1
    Vd3o(3*k-2:3*k,1) = Vg3;
    Sd3(3*k-2:3*k,1) = [Nodos(k+1,6)+1j*Nodos(k+1,7);
                        Nodos(k+1,8)+1j*Nodos(k+1,9);
                        Nodos(k+1,10)+1j*Nodos(k+1,11)];
end
for t = 1:tmax
    for k = 2:NN
        if Nodos(k,12) == 0
            Id3(3*k-5:3*k-3,1) =... 
            -inv(diag(conj(Vd3o(3*k-5:3*k-3,1))))*...
            conj(Sd3(3*k-5:3*k-3,1));
        else
            Id3(3*k-5:3*k-3,1) =...
            -(inv(diag(conj(M*Vd3o(3*k-5:3*k-3,1)))) +...
            inv(diag(conj(M.'*Vd3o(3*k-5:3*k-3,1))))*H)*...
            conj(Sd3(3*k-5:3*k-3,1));
        end
    end
    Vd3 = -IYdd3*(Ydg3*Vg3 - Id3);
    if max(abs(abs(Vd3) - abs(Vd3o))) < epsilon
        Vr3 = [Vg3;Vd3];
        break
    else
        Vd3o = Vd3;
    end
end
%% Impresión de resultados
Ploss = real(Vr3.'*conj(Y3*Vr3))*Sbase;
fprintf('Numero de iteraciones = %i\n',t);
fprintf('\nVoltajes en los nodos por fase(pu)\n');
for n = 1: NN
    fprintf ('\nV%dA=%.4f<%.4f|V%dB=%.4f<%.4f|V%dC=%.4f<%.4f\n',...
    n,abs(Vr3(3*n-2,1)),angle(Vr3(3*n-2,1))*180/pi,...
    n,abs(Vr3(3*n-1,1)),angle(Vr3(3*n-1,1))*180/pi,n,...
    abs(Vr3(3*n,1)),angle(Vr3(3*n,1))*180/pi)
end
fprintf('\nPerdidas de potencia en el sistema(kW)\n');
fprintf('\nPloss = %.4f \n',Ploss);

%% Funciones
function [Zm] = MatrizZ22(C)
if C == 1
    Zm = [1.4923+0.5687i 0.0403+0.2743i 0.0403+0.2864i;
          0.0403+0.2743i 1.4923+0.5687i 0.0403+0.2702i;
          0.0403+0.2864i 0.0403+0.2702i 1.4923+0.5687i];
end
end