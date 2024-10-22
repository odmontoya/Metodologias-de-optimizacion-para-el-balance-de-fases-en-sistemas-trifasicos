%% ALGORITMO DE SENOS Y COSENOS
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Parametros del algoritmo
Ni = 10; tmax = 1000; Nv = 7; ymin = 1; ymax = 6; 

%% Poblacion inicial
p = zeros(Ni,Nv); p = randi([ymin ymax],Ni,Nv); p = [p,zeros(Ni,1)];

%% Algoritmo de senos y cosenos
for t = 0:tmax
    if t == 0
        for k = 1:Ni
            p(k,end) =  N8FP(p(k,1:Nv));
        end
        p = sortrows(p,Nv+1);
        pbest = p(1,:);
    else
        for k = 1:Ni
            r1 = rand; r2 = 1 - t/tmax; r3 = 2*pi*rand; r4 = rand;
            if  r1 < 0.5
                y = round( p(k,1:Nv) + r2*sin(r3)*...
                    abs(r4*pbest(1,1:Nv) - p(k,1:Nv)));
            else
                y = round( p(k,1:Nv) + r2*cos(r3)*...
                    abs(r4*pbest(1,1:Nv) - p(k,1:Nv)));
            end
            y = max(y,ymin);
            y = min(y,ymax);
            y = [y,zeros(1,1)];
            y(1,end) = N8FP(y(1,1:Nv));
            if y(1,end) < p(k,end)
                p(k,:) = y;
            end 
        end
        p = sortrows(p,Nv+1);
        pbest = p(1,:); 
        fprintf('Iteracion = %d Ploss = %.4f\n',t,pbest(1,end));
    end
end