%% MÉTODO DE FLUJO DE POTENCIA POR APROXIMACIONES SUCESIVAS
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.
function [Ploss] = N8FP(Cs)
%% Datos del sistema de prueba
Sbase = 1000; Vbase = 11/sqrt(3); Zbase = ((1000*Vbase)^2)/(1000*Sbase);
Lineas = [1 2 1 1; 2 3 2 1; 2 5 3 1;
          2 7 3 1; 3 4 4 1; 3 8 5 1;
          5 6 6 1];
C1 = exp(1j*-2*pi/3); C2 = exp(1j*2*pi/3);
Nodos = [1 0 1 C1 C2 0   0   0    0  0   0   0 1;
         2 1 1 C1 C2 519 250 259 126 515 250 0 1;
         3 1 1 C1 C2 0	 0	 259 126 486 235 0 1;
         4 1 1 C1 C2 0	 0	 0   0   324 157 0 1;
         5 1 1 C1 C2 0	 0	 0   0   226 109 0 1;
         6 1 1 C1 C2 0	 0	 0   0   145 70  0 1;
         7 1 1 C1 C2 486 235 0   0   0   0   0 1;
         8 1 1 C1 C2 0	 0	 267 129 0   0   0 1];  
Nodos(2:end,13) = Cs'; 
Nodos(:,6:11) = Nodos(:,6:11)/Sbase;       
%% Cálculos iniciales
NN = size(Nodos,1); NL = size(Lineas,1); Y3 = zeros(3*NN,3*NN); 
for i = 1:NL
    Ni = Lineas(i,1); Nj = Lineas(i,2);
    ZL = Matriz8Z(Lineas(i,3))/Zbase;
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
for i = 2:NN
    if Nodos(i,13) == 2 %BCA
        Sd3(3*i-5:3*i-3,1) = [0 1 0; 0 0 1; 1 0 0]*Sd3(3*i-5:3*i-3,1);
    elseif Nodos(i,13) == 3 %CAB
        Sd3(3*i-5:3*i-3,1) = [0 0 1; 1 0 0; 0 1 0]*Sd3(3*i-5:3*i-3,1);
    elseif Nodos(i,13) == 4 %ACB
        Sd3(3*i-5:3*i-3,1) = [1 0 0; 0 0 1; 0 1 0]*Sd3(3*i-5:3*i-3,1);
    elseif Nodos(i,13) == 5 %CBA
        Sd3(3*i-5:3*i-3,1) = [0 0 1; 0 1 0; 1 0 0]*Sd3(3*i-5:3*i-3,1);
    elseif Nodos(i,13) == 6 %BAC
        Sd3(3*i-5:3*i-3,1) = [0 1 0; 1 0 0; 0 0 1]*Sd3(3*i-5:3*i-3,1);
    end
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
end