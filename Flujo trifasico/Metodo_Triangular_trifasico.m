%% MÉTODO DE FLUJO DE POTENCIA POR METODO TRIANGULAR
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Datos del sistema de prueba
clc, clear
Sbase = 1000; Vbase = 13.8; Zbase = ((1000*Vbase)^2)/(1000*Sbase); 
Lineas = [1 2 1; 1 3 1];
C1 = exp(1j*-2*pi/3) ; C2 = exp(1j*2*pi/3) ;
Nodos = [1 0 1 C1 C2 0    0    0    0    0    0    0;
         2 1 1 C1 C2 1600 990  2000 1250 2300 1500 0;
         3 1 1 C1 C2 2400 1500 1300 720  1900 1200 0];        
Nodos(:,6:11) = Nodos(:,6:11)/Sbase;
%% Cálculos iniciales
NN = size(Nodos,1); NL = size(Lineas,1); 
T3 = zeros(3*NL,3*NN); Zp3 = zeros(3*NL,3*NL);
for i = NL:-1:1 
    Ni = Lineas(i,1); Nj = Lineas(i,2);
    T3(3*i-2:3*i,3*Nj-2:3*Nj) = [1 0 0; 0 1 0; 0 0 1];
    for j = NL-1:-1:1
        if Lineas(j,2) == Ni
           T3(3*j-2:3*j,3*Nj-2:3*Nj) = [1 0 0; 0 1 0; 0 0 1];
           Ni = Lineas(j,1); 
        end
    end
    Zp3(3*i-2:3*i,3*i-2:3*i) = MatrizZ26(Lineas(i,3))/Zbase;
end
T3(:,1:3) = [];
Zbb3 = T3.'*(Zp3*T3);
%% Método de flujo de potencia por metodo triangular
tmax = 100; tol = 1e-10;
Vg3 = [Nodos(1,3),Nodos(1,4),Nodos(1,5)].'; 
Vd3o = zeros(3*NN-3,1); unos = ones(size(Vd3o,1),1);
Sd3 = zeros(3*NN-3,1); Id3 = zeros(3*NN-3,1);
H = [0 0 1; 1 0 0; 0 1 0]; M = [1 -1 0; 0 1 -1; -1 0 1];
for k = 1:NN-1
    Vd3o(3*k-2:3*k,1) = Vg3;
    unos(3*k-2:3*k,1) = Vg3;
    Sd3(3*k-2:3*k,1) = [Nodos(k+1,6)+1j*Nodos(k+1,7);
                        Nodos(k+1,8)+1j*Nodos(k+1,9);
                        Nodos(k+1,10)+1j*Nodos(k+1,11)];
end
for t = 1:tmax
    for k = 2:NN
        if Nodos(k,12) == 0 
            Id3(3*k-5:3*k-3,1) =...
            -inv(diag(conj(Vd3o(3*k-5:3*k-3,1))))...
            *conj(Sd3(3*k-5:3*k-3,1));
        else 
            Id3(3*k-5:3*k-3,1) =...
            -(inv(diag(conj(M*Vd3o(3*k-5:3*k-3,1)))) +...
            inv(diag(conj(M.'*Vd3o(3*k-5:3*k-3,1))))*H)*...
            conj(Sd3(3*k-5:3*k-3,1));
        end
    end
    Vd3 = unos + Zbb3*Id3;
    if max(abs(abs(Vd3) - abs(Vd3o))) < tol
        Vr3 = [Vg3;Vd3];
        break
    else
        Vd3o = Vd3;
    end
end
%% Impresión de resultados
J3 = -T3*Id3;
Ploss = real(J3.'*Zp3*conj(J3))*Sbase;
fprintf('Número de iteraciones = %i\n',t);
fprintf('\nVoltajes en los nodos por fase(pu)\n');
for n = 1: NN
    fprintf ('\nV%dA=%.4f<%.4f|V%dB=%.4f<%.4f|V%dC=%.4f<%.4f\n',...
    n,abs(Vr3(3*n-2,1)),angle(Vr3(3*n-2,1))*180/pi,...
    n,abs(Vr3(3*n-1,1)),angle(Vr3(3*n-1,1))*180/pi,n,...
    abs(Vr3(3*n,1)),angle(Vr3(3*n,1))*180/pi)
end
fprintf('\nPérdidas de potencia en el sistema(kW)\n');
fprintf('\nPloss = %.4f \n',Ploss);
%% Funciones
function [Zm] = MatrizZ26(C)
if C == 1
    Zm = [1.5025+0.5088i 0.0287+0.4088i 0.0302+0.3088i;
          0.0287+0.4088i 1.5025+0.5088i 0.0287+0.4088i;
          0.0302+0.3088i 0.0287+0.4088i 1.5025+0.5088i];
end
end