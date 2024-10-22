%% ALGORITMO DE BUSQUEDA POR CUERVOS
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Parametros del algoritmo
Ni = 10; tmax = 1000; Nv = 7; ymin = 1; ymax = 6; Ap = 0.1; fl = 2;

%% Poblacion inicial
p = zeros(Ni,Nv);  
pnew = zeros(Ni,Nv); 
p = randi([ymin ymax],Ni,Nv);
p = [p,zeros(Ni,1)]; 

%% Algoritmo de búsqueda por cuervos
for t = 0:tmax
    if t == 0
        for k = 1:Ni
            p(k,end) =  N8FP(p(k,1:Nv));
        end
        mem = p;
    else
        num = ceil(Ni*rand(1,Ni));
        for i = 1:Ni
            if rand > Ap
                pnew(i,1:Nv) = round(p(i,1:Nv) + fl*rand*...
                               (mem(num(i),1:Nv) - p(i,1:Nv)));
            else
                pnew(i,1:Nv) = randi([ymin ymax],1,Nv);
            end
        end
        pnew = [pnew(:,1:Nv),zeros(Ni,1)];
        for k = 1:Ni
            pnew(k,end) =  N8FP(pnew(k,1:Nv));
        end
        for i=1:Ni 
            if pnew(i,1:Nv) >= ymin & pnew(i,1:Nv) <= ymax
                p(i,1:Nv) = pnew(i,1:Nv); 
                if pnew(i,end) < mem(i,end)
                    mem(i,1:Nv) = pnew(i,1:Nv); 
                    mem(i,end) = pnew(i,end);
                end
            end
        end
        mem = sortrows(mem,Nv+1);
        xbest = mem(1,:);
        fprintf('Iteracion = %d Ploss = %.4f\n',t,xbest(1,end));
    end
end
xbest(1,1:Nv) 