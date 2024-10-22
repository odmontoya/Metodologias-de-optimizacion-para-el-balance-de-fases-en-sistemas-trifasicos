%% ALGORITMO DE AGUJEROS NEGROS
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Parametros del algoritmo
Ni = 10; tmax = 1000; Nv = 7; ymin = 1; ymax = 6; 

%% Poblacion inicial
p = zeros(Ni,Nv); p = randi([ymin ymax],Ni,Nv); p = [p,zeros(Ni,1)]; 

%% Algoritmo de Agujeros Negros
for t = 0:tmax
    if t == 0
        for k = 1:Ni
            p(k,end) =  N8FP(p(k,1:Nv));
        end
        p = sortrows(p,Nv+1);
        BH = p(1,:);
    else
        for k = 2:Ni
            p(k,1:Nv) = round(p(k,1:Nv) + rand*(BH(1,1:Nv) -...
                        p(k,1:Nv))); 
            p(k,1:Nv) = max(p(k,1:Nv),ymin);
            p(k,1:Nv) = min(p(k,1:Nv),ymax);
        end
        for k = 2:Ni
            p(k,end) =  N8FP(p(k,1:Nv));
        end
        p = sortrows(p,Nv+1);
        if p(1,end) < BH(1,end)
            BH = p(1,:);
        end
        R = BH(1,end)/sum(p(:,end));
        for k = 2:Ni
            dist = norm(BH(1,1:Nv) - p(k,1:Nv));
            if dist < R
                p(k,1:Nv) = randi([ymin ymax],1,Nv);
                p(k,end) =  N8FP(p(k,1:Nv));
            end
        end
        p = sortrows(p,Nv+1);
        if p(1,end) < BH(1,end)
            BH = p(1,:);
        end
        fprintf('Iteracion = %d Ploss = %.4f\n',t,BH(1,end));
    end
end
BH(1,1:Nv) 