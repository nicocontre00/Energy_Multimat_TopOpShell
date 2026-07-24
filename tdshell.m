%**************************************************************************
% Topological Derivative for the Compliance  
%**************************************************************************
% DESCRIPTION
% Computes the topological derivative
% 
% INPUT
% mesh:    pdetool mesh struct
% U:       pdetool solution
% matprop: material properties struct
% psi:     level-set function
% params:  topology optimization parameters struct
% pdecoef: pdetool coeficiets struct
%
% OUTPUT
% dtE: compliance T. D. -> dtE = -DT (bulk) and dt = DT (inclusion)
%
% HISTORY
% A.A. Novotny     05/2012: code implementation.
% C.G. Lopes       07/2017: final code updating.
% R.B. dos Santos  07/2017: final code updating.
% S.M. Giusti      08/2018: implementation of RM plate.
%**************************************************************************

function [dt] = tdshell(mesh,U,pdecoef, matprop,signatures,psi,vol,params,comp0,vol0)
penalization = params.penalization;
np = mesh.np; p = mesh.p; t = mesh.t;
h = matprop.h0; % Esto ahora es un vector columna
nu = matprop.nu0;


%% The membrane problem
    % Material properties
    la = matprop.la_andes;
    mu = matprop.mu_andes;
    gamma = matprop.gamma;
    tgamma = pdeintrp(p,t,(psi<0)+gamma*(psi>=0)); 
    
    % Stress
    id=[1;1;0];
    alpha = (la+mu)/mu;
    beta = (la+3*mu)/(la+mu);
    e = e_andes(mesh,U,h,signatures); % strain     
    e(3,:)=1/2.* e(3,:); 
    divu = e(1,:) + e(2,:);  % e(1,:) = exx ; e(2,:) = eyy ; e(3,:) = exy de ANDES
    s = la*id*divu + 2*mu*e; % Nominal stress                
    tgamma3 = [tgamma;tgamma;tgamma];
    s = s.*tgamma3; % Effective stress      
    %Incluiding the body force contribution
    %f = f0;
    %f=0;
    %tU = pdeintrp(p,t,U); 
    %dtf = (f(1,:).*tU(1,:)+f(2,:).*tU(2,:)); 

    % Topological derivative at the bulk phase
    coef0 = 0.5*((1.0-gamma)/(1+gamma*beta)); 
    coef1 = coef0*(1+beta); 
    coefx = (1.0-gamma)/(1.0+gamma*alpha);         
    coef2 = 0.5*coef0*(alpha-beta)*coefx;         
    %    e(1)=exx; e2=exy; e3=eyy; donde e(1,:) = eandes(1,:)
    %                                    e(2,:) = eandes(3,:)
    %                                    e(3,:) = eandes(2,:)
    dte_m = coef1*((s(1,:).*e(1,:)+2*s(3,:).*e(3,:)+s(2,:).*e(2,:))) ...
            + coef2*((s(1,:)+s(2,:)).*(e(1,:)+e(2,:)));% + dtfe ;  

    % Topological derivative at the inclusion
    gamma = 1.0/gamma; 
    coef0 = 0.5*((1-gamma)/(1+gamma*beta));
    coef1 = coef0*(1+beta);
    coefx = (1.0-gamma)/(1.0+gamma*alpha);
    coef2 = 0.5*coef0*(alpha-beta)*coefx;
    dti_m = coef1*((s(1,:).*e(1,:)+2*s(3,:).*e(3,:)+s(2,:).*e(2,:))) ...
        + coef2*((s(1,:)+s(2,:)).*(e(1,:)+e(2,:)));% + dtfi; 
        
%% The Bending problem
    % Material properties
    gamma = matprop.gamma;
    tgamma = pdeintrp(p,t,(psi<0)+gamma*(psi>=0));
    h0 = matprop.h0; % Ahora es un vector columna
    %h3 = h0*h0*h0;
    h3 = h0.*h0.*h0; % Esto da un vector de la misma dimension pero en donde cada elemento esta elevado al cubo
    la=matprop.la_dkt;
    mu=matprop.mu_dkt;

    % Stress
    alpha = (la+mu)/mu;
    beta = mu/(3.0*mu+2.0*la);
    nu = la/(2*mu+la);
    E = 2*mu*(1+nu);
    tE = E*tgamma;
    s = getstress_plate(mesh,U,matprop); % nominal stress
    tgamma3 = [tgamma;tgamma;tgamma]; 
    s = s.*tgamma3; % effective stress

    % Topological derivative at the bulk phase
    aa = s(1,:).*s(1,:) + 2*s(2,:).*s(2,:) + s(3,:).*s(3,:); % s.s
    bb = (s(1,:) + s(3,:)).*(s(1,:) + s(3,:)); % (tr(s))^2
    % coef0 = 6.0*(1.0 - gamma)./(tE*h3);
    coef0 = 12.0*(1.0-gamma)./(tE.*h3'); 
    coef1 = 4.0*(alpha*beta)/(1.0 + beta*gamma);
    coef2 = 1.0/(1.0 + alpha*gamma) - 2.0*(alpha*beta)/(1.0 + beta*gamma);
    dte_b = coef0.*(coef1*aa + coef2*bb);

    % Topological derivative at the inclusion
    gamma = 1.0/gamma;
    % coef0 = 6.0*(1.0 - gamma)./(tE*h3);%% EN EL ORI EL COEF0 ES 12*1-GAMMA...
    coef0 = 12.0*(1.0-gamma)./(tE.*h3'); %% AcA ESTA IGUAL AL COMP ORI
    coef1 = 4.0*(alpha*beta)/(1.0 + beta*gamma);
    coef2 = 1.0/(1.0 + alpha*gamma) - 2.0*(alpha*beta)/(1.0 + beta*gamma);
    dti_b = coef0.*(coef1*aa + coef2*bb);

%% Include the volume constraint sensitivity

    penalty = params.penalty;
	vol0 = params.vol0;
	auglag= params.auglag;
    volfrac = params.volfrac;
	voltarget = vol0*volfrac;

    % topological derivative at the bulk fase. membrane
    dte_m = dte_m / comp0; 
    if penalization == 1 
        dte_m = dte_m - penalty / vol0; 
    elseif penalization == 2
        coef = vol / voltarget;
        dte_m = dte_m - penalty*((coef <= 1)*0.0 - (coef > 1)*2.0*(1.0-coef)) / voltarget;
    elseif penalization == 3
        coef = vol / voltarget;
        dte_m = dte_m - max(0,penalty - auglag*(1.0-coef)) / voltarget;
    end

    % topological derivative at the inclusion. membrane
    dti_m = dti_m / comp0; 
    if penalization == 1 
        dti_m = dti_m + penalty / vol0;                
    elseif penalization == 2
        coef = vol / voltarget;
        dti_m = dti_m + penalty*((coef <= 1)*0.0 - (coef > 1)*2.0*(1.0-coef)) / voltarget;
    elseif penalization == 3
        coef = vol / voltarget;
        dti_m = dti_m + max(0,penalty - auglag*(1.0-coef)) / voltarget;
    end

    % topological derivative at the bulk fase. Bending
    dte_b = dte_b / comp0; 
    if penalization == 1 
        dte_b = dte_b - penalty / vol0; 
    elseif penalization == 2
        coef = vol / voltarget;
        dte_b = dte_b - penalty*((coef <= 1)*0.0 - (coef > 1)*2.0*(1.0-coef)) / voltarget;
    elseif penalization == 3
        coef = vol / voltarget;
        dte_b = dte_b - max(0,penalty - auglag*(1.0-coef)) / voltarget;
    end

    % topological derivative at the inclusion. Bending
    dti_b = dti_b / comp0; 
    if penalization == 1 
        dti_b = dti_b + penalty / vol0;                
    elseif penalization == 2
        coef = vol / voltarget;
        dti_b = dti_b + penalty*((coef <= 1)*0.0 - (coef > 1)*2.0*(1.0-coef)) / voltarget;
    elseif penalization == 3
        coef = vol / voltarget;
        dti_b = dti_b + max(0,penalty - auglag*(1.0-coef)) / voltarget;
    end
	
    %% Suma de ambas derivadas
    % Smoothing of the topological derivative
    pchi = (psi<0); tchi = pdeintrp(p,t,pchi);
    % Compute function g: g = -DT (bulk) and g = DT (inclusion)
    dtE_m = -tchi.*dte_m + (1-tchi).*dti_m; dtE_m = pdeprtni(p,t,dtE_m);
    dtE_b = -tchi.*dte_b + (1-tchi).*dti_b; dtE_b = pdeprtni(p,t,dtE_b);
    dtE_shell = dtE_m + dtE_b;
    %salida para comp shell
    dt=dtE_shell;

end