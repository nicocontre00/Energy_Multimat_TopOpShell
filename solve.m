    %**************************************************************************
% Solve the pdetool FEM linear system  
%**************************************************************************
% DESCRIPTION
% Assemble and solve the pdetool fineite element linear system 
% 
% INPUT
% psi:     level-set function
% mesh:    pdetool mesh struct
% matprop: material properties struct
% pdecoef: pdetool coeficiets struct
% bc:      extra boundry conditions struct
%
% OUTPUT
% [U,F,vol]: pdetool solution, pdetool load and current volume
%
% HISTORY
% J-M.C. Farias     12/2010: code implementation.
% A.A. Novotny      12/2010: code updating.
%**************************************************************************

function [U,F,vol] = solve(psi, mesh, matprop, pdecoef, bc, signatures)
    gamma  = matprop.E./matprop.E(1); % contrast for each material
    p = mesh.p; e = mesh.e; t = mesh.t; area = mesh.area;
    b = pdecoef.b; a = pdecoef.a; f = pdecoef.f; c0 =  pdecoef.c(:,1);
    % element characteristic function: 
    [ ~,tfi ] = charfunc(p,t,psi,pdecoef);
%     tgamma = pdeintrp(p,t,fi*gamma');  %for P1-projection aproach
    tgamma = gamma*tfi; %for mixed formulation aproach
    c = c0 * tgamma; % effective elasticity tensor 
%     [K,F] = assempde(b,p,e,t,c,a,f);
	if size(f,1)==6;
		[K]= assemashell(mesh,pdecoef, matprop, signatures,c);
		F=sparse(mesh.np*mesh.N,1);		
		% Apply Dirichlet boundary condition   
        np = mesh.np;
        eq = bc.pDir(:,1) + (bc.pDir(:,2)-1)*np; neq = size(eq,1);
        K(eq,:) = 0.0; K(:,eq) = 0.0; K(eq,eq) = speye(neq);
        % Apply Neumann boundary condition
        eq = bc.pNeu(:,1) + (bc.pNeu(:,2)-1)*np; 
        F(eq) = F(eq) + bc.pNeu(:,3);   
	else
	    [K,~,F] = assema(p,t,c,a,f);
		[K,F] = pdeupdate(K,F,bc,mesh);
	end
    U = K \ F;  % solve linear system
%     vol  = area*tchi';  %calculate volume assigned to each material including void -- P1 projection aproach
    vol  = area*tfi'; %calculate volume assigned to each material including void -- mixed formulation aproach
end