%% ALGORITMO GENÉTICO DE CHU & BEASLEY
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Parametros del algoritmo
Ni = 10; tmax = 1000; Nv = 7; ymin = 1; ymax = 6; 

%% Poblacion inicial
p = zeros(Ni,Nv); p = randi([ymin ymax],Ni,Nv); 

%% Algoritmo genetico de Chu & Beasle
for k = 1:Ni 
    for j = 1:Ni
        if j ~= k
            while (p(k,1:Nv) == p(j,1:Nv))
                p(j,1:Nv) = randi([ymin ymax],1,Nv);
            end
        end
    end
end
p = [p,zeros(Ni,1)];
for t = 0:tmax
    if t == 0
        for k = 1:Ni
            p(k,end) =  N8FP(p(k,1:Nv));
        end
        p = sortrows(p,Nv+1);
    else
        s = randi([1 Ni],2,1);
        while size(unique(s),1) == 1
            s = randi([1 Ni],2,1);
        end
        P1 = p(s(1),:);
        P2 = p(s(2),:);
        r = randi([1 Nv-1],1,1);
        H1 = [P1(1,1:r), P2(1,r+1:end)];
        H2 = [P2(1,1:r), P1(1,r+1:end)];
        m = randi([1,Nv],1,1);
        H1(1,m) = randi([ymin ymax],1,1);
        H2(1,m) = randi([ymin ymax],1,1);
        ph = [H1;H2];
        for k = 1:2
            ph(k,end) =  N8FP(ph(k,1:Nv));
        end
        ph = sortrows(ph,Nv+1);
        if ph(1,end) < p(end,end)
            ban = 0;
            for k = 1:Ni
                if p(k,1:Nv) == ph(1,1:Nv)
                    ban = 1;
                    break;
                end
            end
            if ban == 0
                p(end,:) = ph(1,:);
            end
        end
        p = sortrows(p,Nv+1); 
        fprintf('Iteracion = %d Ploss = %.4f\n',t,p(1,end));
    end
end
p(1,1:Nv)