%% ALGORITMO DE BUSQUEDA POR VORTICES
%% Libro: Metodologías de optimización para el balance de fases en sistemas trifásicos
%% Autores: L. S. Avellaneda-Gómez, B. Cortés-Caicedo, O. D. Montoya-Giraldo
%% Editorial Universidad Distrital Francisco José de Caldas, 2024.

%% Parametros del algoritmo
Ni = 10; tmax = 1000; Nv = 7; ymin = 1; ymax = 6; w = 0.1;

%% Poblacion inicial
mu = ((ymin + ymax)/2)*ones(1,Nv); mu = round(mu);
sigma = (ymax - ymin)/2; sigma = round(sigma); r = sigma;
Pbest0 = inf;

%% Algoritmo de búsqueda por vortices
for t = 1:tmax
    P = r*randn(Ni,Nv);
    P = round(P + mu);
    P = max(P,ymin+0.5);
    P = min(P,ymax-0.5);
    P = [P,zeros(Ni,1)];
    for k = 1:Ni
        P(k,end) =  N8FP(P(k,1:Nv));
    end
    P = sortrows(P,Nv+1);
    Pbest = P(1,:);
    if Pbest(1,end) < Pbest0
        Pbest0 = Pbest(1,end);
        mubest = Pbest(1,1:Nv);
    end
    mu = mubest;
    at = 1 - t/tmax;
    r = sigma*(1/w)*gammaincinv(w,at);
    fprintf('Iteracion = %d Ploss = %.4f\n',t,Pbest0);
end
Pbest(1,1:Nv)