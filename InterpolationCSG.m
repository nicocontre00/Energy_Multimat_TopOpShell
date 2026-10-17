function [psi] = InterpolationCSG (psi,mesh,msh)

%% CORDENADAS BARICENTRICAS

t_old = mesh.t; p_old = mesh.p; rg_old = mesh.RG';

id1 = t_old(1,:); id2 = t_old(2,:); id3 = t_old(3,:);

x1 = p_old(1,id1); x2 = p_old(1,id2); x3 = p_old(1,id3);
y1 = p_old(2,id1); y2 = p_old(2,id2); y3 = p_old(2,id3);
z1 = p_old(3,id1); z2 = p_old(3,id2); z3 = p_old(3,id3);

x_g = (x1+x2+x3)/3; y_g = (y1+y2+y3)/3; z_g = (z1+z2+z3)/3;

G = [x_g; y_g; z_g]; 

p_new = msh.POS'; t_new = msh.TRIANGLES'; t_new(4, :) = []; rg_new = msh.RG';

%% PREPARACION

G_rg = [G;rg_old]; 
t_rg = [t_new;rg_new]; 
num_rg = size(unique(rg_old),2);
num_nodes = size(p_new,2); 

%% ASIGNO VALOR DE REGION GEOMETRICA A LOS NODOS NUEVOS (INCLUSIVE LOS DE FRONTERA)

nodes = []; % Creo la variable vacia

for j = 1:num_nodes % Repito el bucle como tantos nodos tenga
	
	[~, columnas] = find(t_rg(1:3, :) == j); % Extraigo los ID de los elementos (columnas) que contienen al nodo j

	rg_j = unique(t_rg(4, columnas)); % Tomo solo los elementos de la cuarta fila que se encuentren en las columnas encontradas y defino los valores unicos en un vector columna
	
	if length(rg_j) > 1 % Si tiene mas de 1 region unica (Que no se repite) ingreso
		
		p_rgj = p_new(:,j); % Armo un vector columna que contiene las Coord. del nodo j (La columna j y las 3 filas)
		
		G_rgj = G_rg(:,ismember(G_rg(4,:), rg_j)); % Armo una matriz que contiene las cordenadas baricentricas y su ID de region geometrica de las regiones del nodo j , me da todas las filas en donde los valores de rg_j que coinciden con la cuarta fila de G_rg 
		
		[~,id] = min((p_rgj(1) - G_rgj(1,:)).^2 + (p_rgj(2) - G_rgj(2,:)).^2 + (p_rgj(3) - G_rgj(3,:)).^2); % Me devuelve el indice del elemento cuya distancia es la menor			
        
		rgj = G_rgj(4,id); % Le asigno el valor de la region geometrica en cuyo triangulo sea el mas cercano, busco en la cuarta fila el elemento de la columna id
	
		nodes = [nodes, [j; rgj]]; % Guardo en una matriz el valor del ID del nodo y su valor correspondiente ID de region geometrica, en cada fila un nodo y su ID
	
	else
	
		nodes = [nodes, [j; rg_j]]; % Sino guardo solo el elemento y su valor correspondiente de region geometrica, en cada columna un nodo y su ID
	end
	
end

p_rg = [p_new;nodes(2,:)]; % Armo una matriz que tiene las cordenadas nuevas y su correspondiente valor de region geometrica

%% INTERPOLACION PARA ELEMENTOS

psi_new = zeros(num_nodes,size(psi,2));

for e = 0:num_rg-1 % Repito el bucle como tantas regiones geometricas tenga
	
	index = find(G_rg(4,:) == e); % Obtengo los indices de la matriz de cordenadas baricentricas que cumple pertenecen a la region e, indices de columnas
	
	G_rge = G_rg(1:3,index); % Armo una matriz que solo tiene las Coord. Baricentricas de los elementos de la region e, cada columna tiene un elemento
	
    index2 = find(p_rg(4,:) == e); % Obtengo los indices de la matriz de cordenadas que cumple pertenecen a la region e, indices de columnas
	
	p_rge = p_rg(1:3,index2); % Armo una matriz que solo tiene las Coord. de los nodos de la region e, un nodo por columna
	
	num = size(p_rge,2); % numero de nodos de la region e, cuento las columnas
	
	for i = 1:num % Inicio el bucle para la interpolacion, tomando cada elemento de la matriz p_rge de cordenadas de la region e
        
		x_new = p_rge(1,i); y_new = p_rge(2,i); z_new = p_rge(3,i); % Extraigo las cordenadas del nodo i
		
		[~,id] = min((x_new - G_rge(1,:)).^2 + (y_new - G_rge(2,:)).^2 + (z_new - G_rge(3,:)).^2); %Me devuelve el indice del elemento cuya distancia es la menor
		
		v1 = p_old(:,t_old(1,index(id))); v2 = p_old(:,t_old(2,index(id))); v3 = p_old(:,t_old(3,index(id))); %Obtengo las cordenadas de los vertices del elemento cuya distancia es menor
		
		P = [x_new; y_new; z_new]; % Vector de Cordenadas del nodo new

		n = cross(v2-v1,v3-v1); % Vector normal

		d = dot(n,P-v1)/norm(n); % Distancia 

		pr = P - d*n/norm(n); % Proyeccion 

		% ahora necesito calcular las funciones de forma para interpolar

		A = 0.5*norm(cross(v2-v1,v3-v1));
		A1 = 0.5*norm(cross(pr-v2,v3-v2));
		A2 = 0.5*norm(cross(pr-v1,v3-v1));
		A3 = 0.5*norm(cross(pr-v1,v2-v1));

		%Armo el vector de funciones de forma
		N = [A1/A; A2/A; A3/A];
		
		psi_v = psi([t_old(1,index(id)); t_old(2,index(id)); t_old(3,index(id))],:);
		psi_new(index2(i),:) = N.' * psi_v;
	end
	
end

psi = psi_new;

end