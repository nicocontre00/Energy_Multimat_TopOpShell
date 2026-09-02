function [mesh, pdecoef, matprop, params, bc, psi, signatures] = shcog_2m(mesh, psi)
%% Mesh selection

meshes = {'shcog_msh2.m','shcog_msh3.m','shcog_msh4.m'};

if isempty(mesh)
    index = 1;
    run(meshes{index});	
    t = msh.TRIANGLES';
else
	index = mesh.index;
	index = index + 1;
    run(meshes{index});    
    cd ..  
    cd ..
	    [psi]=Interpolation(psi,mesh,msh);
    cd('examples\shell')		
end

cd ..\..
%% Mesh struct
[mesh]=meshstruct2(msh);        

mesh.index = index;

% Rotation of elements   
[mesh]=rotation(mesh);

%% Material properties

gamma = 1.0e-4;
nu = 0.3; 
h = 1.1;  
E0 = 210.0E+3;
E = E0*[1 1e-4];

matprop.nu0 = nu;
matprop.h0 = h;
matprop.E = E;
matprop.E0 = E;
    
% DKT (Discrete Kirchhoff Triangle)
la = nu.*E/((1+nu)*(1-2*nu)); 
mu = E./(2*(1+nu)); % plane strain
la = 2.*mu.*la./(la+2.*mu); % plane stress

matprop.la_dkt = la; 
matprop.mu_dkt = mu;
    
% ANDES (Assumed natural deviatoric strain)
matprop.la_andes = ((nu.*E.^2)./((1+nu).^2).*(1-2.*nu))./(((nu.*E)./(1+nu).*(1-2.*nu))+(E./(1+nu)));
matprop.mu_andes = E./(2.*(1+nu));
    
% pde coefficients
[pdecoef] = cprop(matprop);    
    
% ANDES signatures   
[signatures]=andes_signature('OPT',nu);
         
%% parameters

% line-search procedure
params.kmin = 1.0E-3; 
% stop criterion
params.stop = 5.0*pi/180; 

% method of penalization
% (1) linear 
% (2) exact quadratic
% (3) augmented-lagrangian
params.penalization = 3;

% penalty parameter
params.penalty = [0.0];   % method 1,2 and 3	% La penalizacion en este caso es solo un valor inicial
params.volfrac = [0.6];   % method 2 and 3		% Porcentaje del volumen final que buscamos (En este ejemplo es 20%)
params.voleps  = [0.01];  % method 2 and 3		% Tolerancia que se puede admitir en el volumen final 0.01
params.auglag  = [0.7];   % method 3			% Coeficiente de actualizacion del multiplicador de Lagrange (penalty aumenta mas o menos rapido) (Profe recomendo un 0.1)
params.epsilon = 0.01;  % method 3 

%% boundary conditions
bc.pNeu = []; bc.pDir = []; 
x = mesh.p(1,:); y = mesh.p(2,:); z = mesh.p(3,:);
eps = 0.01*sqrt(min(mesh.A));

% Dirichlet
cx = 0; 
node1 = find(x < cx + eps); node1 = node1';
cy = 0;
node2 = find(y < cx + eps); node2 = node2';
cx = 20;     
node3 = find(x > cx - eps); node3 = node3';
cy = 20;     
node4 = find(y > cx - eps); node4 = node4';
node = [node1;node2;node3;node4]; npb = size(node,1);
node = [node;node;node;node;node;node];
dof = [ones(npb,1);2*ones(npb,1);3*ones(npb,1);4*ones(npb,1);5*ones(npb,1);6*ones(npb,1)];
val = [zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1);zeros(npb,1)];
bc.pDir = [node, dof, val];

% Neumman
cx = 10; cy = 10; cz = 10;
[~,node4] = min((z - cz).^2 + (x - cx).^2 + (y-cy).^2);
bc.pNeu = [node4, 3,-10]; 

%% initial guess
if index == 1
    psi = ones(size(mesh.p,2),length(matprop.E)-1); psi(:,1) = -1;	
end

mesh.ghold = [];
pdecoef.f = zeros(6,1);
pdecoef.b = 0;
pdecoef.a = 0;

cd('examples\shell') 
end