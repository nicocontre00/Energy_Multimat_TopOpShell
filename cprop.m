function [pdecoef] = cprop(matprop)
	E = matprop.E0; h = matprop.h0; nu = matprop.nu0;
	mu = matprop.mu_andes; la = matprop.la_andes;
	
	c = zeros(16,length(E));
	
	%% DKT coefficients
	
	c(1,:) = E;
	c(2,:) = (E.*h.^3)./(12.*(1-nu.^2));
	c(3,:) = c(2,1)*nu;

	%% ANDES coefficients
	
	c(4,:) = 2*mu + la;       
	c(5,:) = 0;               
	c(6,:) = mu;              
	c(7,:) = 0;               
	c(8,:) = mu;              
	c(9,:) = la;              
	c(10,:) = 0;              
	c(11,:) = mu;                
	c(12,:) = 0;              
	c(13,:) = 2*mu + la;       

	pdecoef.c = c;
	
%	pdecoef.c_dkt = c(2:3,:);
%	pdecoef.c_andes = c(4:13,:);
end