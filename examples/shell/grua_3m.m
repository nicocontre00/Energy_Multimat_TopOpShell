function [mesh, pdecoef, matprop, params, bc, psi, signatures] = grua_3m

%function [mesh, params, psi0, bc, signatures, pdecoef, matprop] = gruab(mesh,psi,params)
    %%cargar malla
   
% mallas = {'gruab_mallag.m','gruab_mallam.m','gruab_mallaf.m','gruab_mallaf02.m','gruab_mallaf03.m'};

%   if isempty(mesh)
	
%	    index = 1;
%        run(mallas{index});
		
%        t = msh.TRIANGLES';
		
%        psi0=-ones(length(unique(t(1:end-1,:))),1);

%   else
        
%		index = mesh.index;
%		index = index + 1;

%        run(mallas{index});
        
%        cd ..       
%		    [psi]=InterpolationCSG(psi,mesh,msh);
%        cd('examples')		

%        psi0 = psi;
%    end

    %% Mesh struct

	run('gruab_mallag.m');
	
        cd ..  
		cd ..
             [mesh]=meshstruct2(msh);
        cd('examples')        

% 	    mesh.index = index;
        mesh.RG = msh.RG;

    %rotaciones de elementos

        cd ..     
             [mesh]=rotation(mesh);
        cd('examples') 

    %% Material properties
  
%   E = 200.0E+9; 
    gamma = 1.0e-4;
    nu = 0.3; 
    h = 0.012;  


% 	idx_h = 2; % 1 para h cte y 2 para h variable
% 
%     if idx_h == 1
%     	h = 0.012;
%         h = h + (zeros (size(mesh.t,2),1));
%     else
%         rg = msh.RG; rg = rg + 1; % Defino el vector de regiones geometricas
%         rg_h = msh.H; rg_h(:,1) = rg_h(:,1) + 1; % Defino la matriz de espesores por region geometrica
%         h = rg_h(rg, 2); %
%     end

%     matprop.E0 = E;
%     matprop.gamma = gamma;
    matprop.nu0 = nu;
 	matprop.h0 = h;
%     matprop.E1 = gamma*E;
    
    E0 = 210.0E+3;
    E = E0*[1 0.5 1e-4];
    matprop.E = E;
    matprop.E0 = E;
    
    
    %DKT (Discrete Kirchhoff Triangle)
    
    la = nu.*E/((1+nu)*(1-2*nu)); 
    mu = E./(2*(1+nu)); % plane strain
    la = 2.*mu.*la./(la+2.*mu); % plane stress
    
    matprop.la_dkt = la; 
    matprop.mu_dkt = mu;
    
    %%ANDES 
	
    matprop.la_andes = ((nu.*E.^2)./((1+nu).^2).*(1-2.*nu))./(((nu.*E)./(1+nu).*(1-2.*nu))+(E./(1+nu)));
    matprop.mu_andes = E./(2.*(1+nu));
    
    %% Coeficiente c pde

    cd .. 
        [pdecoef] = cprop(matprop);    
    cd('examples')    

    %%Firma andes   

    cd ..
        [signatures]=andes_signature('OPT',nu);
    cd('examples')         
   
    %% parameters
        % line-search procedure
        params.kmin = 1.0E-3; % Paso minimo de busqueda
        % stop criterion
        params.stop = 1.0*pi/180; % Valor del angulo en el cual converge el programa

    % method of penalization
    % (1) linear 
    % (2) exact quadratic
    % (3) augmented-lagrangian
    params.penalization = 3;

    % penalty parameter
%     params.penalty = 0;     % method 1, 2 and 3		% La penalizacion en este caso es solo un valor inicial
%     params.volfrac = 0.95;   % method 2 and 3		% Porcentaje del volumen final que buscamos (En este ejemplo es 20%)
%     params.voleps  = 0.01;  % method 2 and 3		% Tolerancia que se puede admitir en el volumen final 0.01
%     params.auglag  = 150;   % method 3				% Coeficiente de actualizacion del multiplicador de Lagrange (penalty aumenta mas o menos rapido) (Profe recomendo un 0.1)
%     params.epsilon = 0.01;  % method 3 

    params.penalty = [0.0 0.0];   % method 1,2 and 3
    params.volfrac = [0.45 0.45];   % method 2 and 3
    params.voleps  = [0.01 0.01];  % method 2 and 3
    params.auglag  = [50 50];   % method 3
    params.epsilon = 0.01;  % method 3

    %% boundary conditions
    
    bc.pNeu = []; bc.pDir = []; 
    x = mesh.p(1,:); y = mesh.p(2,:); z = mesh.p(3,:);
    eps = 0.01*sqrt(min(mesh.A));

%         cx = 0;     
%         node = find(x < cx + eps); node = node';
%         npb = size(node,1); node = [node;node;node;node;node;node];
%         dof = [ones(npb,1);2*ones(npb,1);3*ones(npb,1);4*ones(npb,1);5*ones(npb,1);6*ones(npb,1)];
%         val = [zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1)];
%         bc.pDir = [node, dof, val];

        cx = 0;     
        node1 = find(x < cx + eps); node1 = node1';
        npb1 = size(node1,1); node1 = [node1;node1;node1;node1;node1];
        dof1 = [ones(npb1,1);2*ones(npb1,1);3*ones(npb1,1);4*ones(npb1,1);5*ones(npb1,1)]; % 5 DOFs
        val1 = [zeros(npb1,1);zeros(npb1,1);zeros(npb1,1);zeros(npb1,1);zeros(npb1,1)];
        bc_1 = [node1, dof1, val1];

        cz = -1.369;     
        node2 = find(z < cz + eps); node2 = node2';
        npb2 = size(node2,1);
        dof2 = [ones(npb2,1)]; % en x
        val2 = [zeros(npb2,1)];
        bc_2 = [node2, dof2, val2];  
        
        cx = 0.848; cy = 0; cz= -1.369;
        [~,node3] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
        
        bc_3 = [node3, 2, 0];
        
        bc.pDir = [ bc_1;
		            bc_2;
					bc_3 ];					
	%% cond de neumman

    cx = 4.15; cy = 0; cz = 0.088;
	
%  	eps1 = 1.01 * min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
%   node1 = min(find((z-cz).^2 + (x-cx).^2 + (y-cy).^2 <= eps1));
    [~,node4] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);

    bc.pNeu = [node4, 3,-15000]; 

    %% initial guess
    psi = ones(size(mesh.p,2),length(matprop.E)-1); psi(:,1) = -1; psi(:,2) = 1;

   mesh.ghold = [];
   pdecoef.f = zeros(6,1);
   pdecoef.b = 0;
   pdecoef.a = 0;
end