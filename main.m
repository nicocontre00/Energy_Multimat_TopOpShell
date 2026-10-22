%**************************************************************************
% Compliance Topology Optimization with Volume Constraint in Plane Stress,
% Steady-State Heat Conduction and Kirchhoff Plates considering multiple
% materials.
%**************************************************************************
%
% DESCRIPTION
% Computes the topological derivative and use it together with a level-set 
% domain representation method in the multi-material topology optimization
% design context.
%
% HISTORY
% S. Amstutz       06/2009: code implementation.
% A.A. Novotny     06/2009: code updating.
% D.E. Campeão     12/2010: code updating.
% J-M.C. Farias
% A.A. Novotny
% A.A. Romero      11/2018: multimaterial extension.
%**************************************************************************
% MAIN REFERENCE
% Onco, A. R., & Giusti, S. M. (2020). 
% A robust topological derivative-based multi-material optimization
% approach: Optimality condition and computational algorithm. Computer
% Methods in Applied Mechanics and Engineering, 366, 113044.
%**************************************************************************

clear all; close all; clc

mesh = []; psi=[];

% Load problem data
% cd('examples\elas')
%     [mesh, pdecoef, matprop, params, bc, psi] = bridge;
%     [mesh, pdecoef, matprop, params, bc, psi] = viga1x2_4m;
%     [mesh, pdecoef, matprop, params, bc, psi] = viga1x2_2m;
%     [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_2m;
%     [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_3m;
%     [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_4m;
%     [mesh, pdecoef, matprop, params, bc, psi] = wall_2m;
% cd('examples\heat')
%     [mesh, pdecoef, matprop, params, bc, psi] = engine;
%     [mesh, pdecoef, matprop, params, bc, psi] = star_2m;
%     [mesh, pdecoef, matprop, params, bc, psi] = star_3m;
%     [mesh, pdecoef, matprop, params, bc, psi] = star_4m;
% cd('examples\plate')
%     [mesh, pdecoef, matprop, params, bc, psi] = cog_3m;
%     [mesh, pdecoef, matprop, params, bc, psi] = cog_4m;
%     [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_2m_mod;
%     [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_2m;
%     [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_3m;
%     [mesh, pdecoef, matprop, params, bc, psi] = symmetry_cog_4m;
%     [mesh, pdecoef, matprop, params, bc, psi] = momento_2m;
%     [mesh, pdecoef, matprop, params, bc, psi] = claDKT2;
%     [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_2m;
%     [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_3m;
%     [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_4m;
cd('examples\shell')
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = grua_2m;
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = grua_3m;
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shviga2x1_2m(mesh, psi);
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shviga2x1_3m(mesh, psi);
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shcog_2m(mesh, psi);
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shcog_3m(mesh, psi);
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = wing_2m(mesh, psi);
%     [mesh, pdecoef, matprop, params, bc, psi, signatures] = wing_3m(mesh, psi);
    [mesh, pdecoef, matprop, params, bc, psi, signatures] = wing_4m(mesh, psi);
cd ..\..
filename = 'wall_2m';
%Configure default settings for plot functions
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 
set(groot,'defaultAxesFontSize',16)

% mesh and geometry  parameter
p = mesh.p; e = mesh.e; t = mesh.t; g = mesh.g;

% topology optimization parameters
stop = params.stop;  kmin = params.kmin; penalization = params.penalization; 
volfrac = params.volfrac; auglag = params.auglag; voleps = params.voleps;

% freeze nodes
phold = [];
f = pdecoef.f;
if isfield(mesh, 'RG')
rg = mesh.RG';
ghold = mesh.ghold';
    for i = 1:size(mesh.ghold,1)
        phold = cat(1,phold,unique(t(1:3,rg == ghold(i))));
    end
else
    for i = 1:size(mesh.ghold,1)
        phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i))));
    end
end
pfree = setdiff(1:size(p,2),phold); % free nodes


f = pdecoef.f;
if size(f,1)==6;
    figure(1); clf; set(1,'WindowStyle','docked'); 
    trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),'FaceColor',[0.5 0.5 0.5]); 
    axis equal
else
    figure(10); clf;
    pdeplot(p,e,t,'xydata',t(4,:),'xystyle','flat','colormap','gray',...
                  'xygrid','off','colorbar','off','title','Geometry'); 
    axis image; axis off;

    figure(1); clf;
    set(1,'WindowStyle','docked');
    pdemesh(p,e,t); axis image; axis off;
end

% shape functional and volume associated to the hold-all domain
psi_hold_all = ones(length(p),length(matprop.E)-1); psi_hold_all(:,1) = -1;

if size(f,1)==6;
else
signatures = 0;
end

[U,F,vol_hold_all] = solve(psi_hold_all, mesh, matprop, pdecoef, bc, signatures); 
if size(f,1)==6;
U = U.U;
end
params.energy0 = 0.5 * dot(F,U); % initial energy 
max_vol = vol_hold_all(1); %volinit = volinit.*gamma; %calculate maximum vol for each fase
params.max_vol = max_vol; % store maximum volume for each fase
voltarget = max_vol.*volfrac;

% stop criterion over the volume constraint 
volstop = voleps.*max_vol;
if penalization == 1
    volstop = max_vol;
end

[~,unitM,~] = assema(p,t,0,1,0); % mass matrix of unity density
psi = psi./normL2( unitM,psi ); % level-set function nomalization 

if size(f,1)==6;
else
% plot hold-all domain
figure(2); clf;
set(2,'WindowStyle','docked');
pdeplot(p,e,t,'xydata',ones(size(p,2),1),'xystyle','flat','colormap','gray',...
              'xygrid','off','colorbar','off'); axis image; axis off; 
end
k = 1; iter = 0;
gsf = []; gEpot = []; gth = []; gvol = [];
option = 'null';

while not(strcmp(option,'s'))
    %-----------Initial values-----------%
    [U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc, signatures);
    [sf,Epot] = shfunc(F,U,volume,params,pdecoef);

    % compute function g: g = -DT (bulk) and g = DT (inclusion)
    dt = topder(mesh,U,volume,matprop,psi,params,pdecoef,signatures);
    dt = dt/normL2(unitM,dt);  % g function normalization

    dt(phold,:) = psi(phold,:); % freeze dt function

    for i = 1 : length(matprop.E)-1
        cosin(i) = psi(:,i)'*unitM*dt(:,i);
    end
    cosin = max(min(sum(cosin),1.0),-1.0);
    theta = max(real(acos(cosin)),1.0e-4);

    difvol = volume(1:end-1)-voltarget;
    ic = (abs(difvol) > volstop); %index control
    while and( or(any(ic),theta > stop) , k/2 > kmin)

        iter = iter + 1;
        [U,F,volume] = solve(psi, mesh, matprop, pdecoef, bc, signatures);
        [sf,Epot] = shfunc(F,U,volume,params,pdecoef);

        % compute function g: g = -DT (bulk) and g = DT (inclusion)
        dt = topder(mesh,U,volume,matprop,psi,params,pdecoef,signatures);
        dt = dt./normL2(unitM,dt);  % g function normalization

        dt(phold,:) = psi(phold,:); % freeze dt function
 
        for i = 1 : length(matprop.E)-1
            cosin(i) = psi(:,i)'*unitM*dt(:,i);
        end
        
        cosin = max(min(sum(cosin),1.0),-1.0);
        theta = max(real(acos(cosin)),1.0e-4);

        % performs a line-search
        sfold = sf; psiold = psi; sf = sf + 1; k = min(1-(1e-5),1.5*k); % trick
%
        while and((sf > sfold) , k>kmin)
%           (sf - sfold > eps(sf))
            % update level-set function
            psi(pfree,:)  = (sin((1-k)*theta)*psiold(pfree,:)...
                +  sin(k*theta)*dt(pfree,:))./sin(theta);
            [U,F,volume] = solve(psi,mesh,matprop,pdecoef,bc,signatures);
            [sf,Epot] = shfunc(F,U,volume,params,pdecoef);
            k = k / 2;
        end
        psi = psi/normL2(unitM,psi);
        k = k * 2;

        [fi,~] = charfunc( p,t,psi,pdecoef); %calculate characteristic function for plot

        gsf  = [gsf,sf];
        gEpot = [gEpot,Epot];
        gth  = [gth,theta];
        gvol = [gvol,volume(1:end-1)'*100/max_vol];

        disp(['iter    = ', num2str(iter)]);
        disp(['volume  = ', num2str(volume(1:end-1)),' => ', sprintf(repmat(' %3.2f%%',1,length(volume)-1),volume(1:end-1)*100/max_vol)]);
        disp(['sf      = ', num2str(sf)]);
        disp(['k       = ', num2str(k)]);
        disp(['theta   = ', num2str(theta*180/pi)]);
        disp(['penalty = ', num2str(params.penalty)]);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

        facecolor = [186 212 244]/255;
        linecolor = [0 62 102]/255;

        if size(f,1)==6;
        colors = [0 0 0;1 0 0; 0 0 1;0 1 0];
        colors = colors(1:size(fi,2),:);
        colors(end,:) = [1 1 1];
        coef = 1:size(fi,2);
        fi_shell = fi*coef';
        figure(3); clf; set(3,'WindowStyle','docked');
        trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),fi_shell,'EdgeColor','none','FaceColor','interp')
        colormap(colors);
        caxis([1 size(fi,2)])
        axis equal
        if isfield(mesh,'RG')
            rg  = mesh.RG(:);
            rib = mesh.rib(:);
            id_rib = ismember(rg,rib);
            nmat = size(fi,2);
            nelem = size(t,2);
            fi_elem = zeros(nelem,nmat);
            for im = 1:nmat
                aux = pdeintrp(p,t,fi(:,im));
                fi_elem(:,im) = aux(:);
            end
            tri_rib = t(1:3,id_rib)';    
            fi_rib  = fi_elem(id_rib,:);  
            fi_solid = fi_rib(:,1:end-1);
            solid_sum = sum(fi_solid,2);
            solid_sum = max(0,min(1,solid_sum));
            tol_vis = 1e-8;
            visible = solid_sum > tol_vis;
            tri_plot   = tri_rib(visible,:);
            fi_plot    = fi_solid(visible,:);
            solid_plot = solid_sum(visible);
            fi_norm = fi_plot ./ solid_plot;
            solidColors = colors(1:end-1,:);
            RGB_solid = fi_norm * solidColors;
            RGB = RGB_solid .* solid_plot + (1 - solid_plot);
            RGB = max(0,min(1,RGB));
            figure(7); clf;
            set(7,'WindowStyle','docked');
        
            patch('Faces',tri_plot, ...
                  'Vertices',p', ...
                  'FaceVertexCData',RGB, ...
                  'FaceColor','flat', ...
                  'EdgeColor','none');
            axis equal
            axis tight
            view(3)
            title('Distribucion multimaterial - Costillas');
        
            drawnow;
        end

        else
        if iter == 1 
            figTopology = figure(3); 
            set(3,'WindowStyle','docked');
            axTopology = axes('Parent', figTopology);
        end
        multimat_plot( p,e,t,fi,axTopology);
        drawnow
        end

%       FRAMES(iter) = getframe(gcf);
        figure(4); clf;hold on, plot(gsf,'k-','LineWidth',1.2); %title('Shape Function');
        ymax = max(max(gsf),max(gEpot));
        yyaxis left; ylabel('Shape Function'), ylim([0 ymax])%ylim('auto') 
        yyaxis right; plot(gEpot,'-','LineWidth',1.2); ylabel('Epot'),grid minor, xlim([0 iter]),ylim([0 ymax]) %ylim('auto') 
        set(4,'WindowStyle','docked');
        figure(5); clf; plot(gth*180/pi,'-','MarkerFaceColor',facecolor,'Color',linecolor,'LineWidth',1.2); title('Theta Angle');
        grid minor , xlim([0 iter]),ylim([min(gth*180/pi)-5 max(gth*180/pi)+10])
        set(5,'WindowStyle','docked');   
        drawnow
        figure(6); clf;
        for index = 1:size(gvol,1)
            plot(gvol(index,:),'-','LineWidth',1.2), hold on
            leyenda{index} = sprintf('$V_{f%d}$',index);
        end
        title('Volume'); legend(leyenda),grid minor, xlim([0 iter]), ylim([-5 105])    
        set(6,'WindowStyle','docked');

        %Update augmented lagrangian parameters
        difvol = volume(1:end-1)-voltarget;
        ic = (abs(difvol) > volstop); %index control

        if any(ic==1) 
            if penalization == 2 % increase the penalization parameter
                params.penalty = 2.0 * params.penalty;
                k = 1;
            elseif penalization == 3 % increase the lagrangian multiplier
                coef = volume(1:end-1)./ voltarget; coef = coef(ic); tau = auglag;
                penalty = params.penalty(ic);
                penalty = penalty + (tau(ic)./auglag(ic)) .* (max(0,penalty - auglag(ic).*(1.0-coef))-penalty);
                params.penalty(ic) = penalty;
                k = 1;
            end
        end
    end

    option = 'null';
    while and(not(strcmp(option,'r')), and(not(strcmp(option,'s')), not(strcmp(option,'c'))))
        option = input('\n -> type "r" to remesh or "s" to stop : ', 's');
    end
    if (option == 'r')
        
        if size(f,1)==6;

            cd('examples\shell')
            %     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shviga2x1_2m(mesh,psi);
            %     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shviga2x1_3m(mesh, psi);
            %     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shcog_2m(mesh, psi);
            %     [mesh, pdecoef, matprop, params, bc, psi, signatures] = shcog_3m(mesh, psi);
            %     [mesh, pdecoef, matprop, params, bc, psi, signatures] = wing_2m(mesh, psi);
%                   [mesh, pdecoef, matprop, params, bc, psi, signatures] = wing_3m(mesh, psi);
                [mesh, pdecoef, matprop, params, bc, psi, signatures] = wing_4m(mesh, psi);
            cd ..\..

            % update data            
            % mesh and geometry  parameters
            p = mesh.p; t = mesh.t;
            area = mesh.A; np = mesh.np; vol0 = sum(area);
            
            % freeze nodes
%             phold = [];
%             for i = 1:size(mesh.ghold,1)
%                 phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i))));
%             end
%             pfree = setdiff(1:size(p,2),phold); % free nodes

phold = [];

rg = mesh.RG(:)';
ghold = mesh.ghold(:)';

for i = 1:length(ghold)
    phold = [phold; unique(t(1:3,rg == ghold(i)))];
end

phold = unique(phold);

pfree = setdiff(1:size(p,2),phold);

            figure(1); clf; set(1,'WindowStyle','docked'); 
            trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),'FaceColor',[0.5 0.5 0.5]); 
            axis equal

            % topology optimization parameters
	
            [~,unitM,~] = assema(p,t,0,1,0); % mass matrix of unity density
            psi = psi./normL2( unitM,psi ); % level-set function nomalization  

            % shape functional and volume associated to the hold-all domain
            psi_hold_all = ones(length(p),length(matprop.E)-1); psi_hold_all(:,1) = -1;
            [U,F,vol_hold_all] = solve(psi_hold_all, mesh, matprop, pdecoef, bc, signatures);
            if size(f,1)==6
                U = U.U;
            end
            params.energy0 = 0.5 * dot(F,U); % initial energy
            
            max_vol = vol_hold_all(1);
            params.volinit = max_vol; % store maximum volume for each fase
            voltarget = max_vol.*volfrac;
            params.max_vol = max_vol;

            k = 1;
            
            option = 'null';

            % stop criterion over the volume constraint
            volstop = voleps.*max_vol;
            if penalization == 1
                volstop = max_vol;
            end
        else

            for i = 1:(strcmp(mesh.remesh,'longest')*2 + strcmp(mesh.remesh,'regular'))
                [p,e,t,psi] = refinemesh(g,p,e,t,psi,mesh.remesh);
                mesh.p = p; mesh.e = e; mesh.t = t;
            end
            
            % update mesh parameters
            figure(1); clf;
            set(1,'WindowStyle','docked');
            pdemesh(p,e,t); axis image; axis off;
            [~,unitM,~] = assema(p,t,0,1,0);
            area = pdetrg(p,t); mesh.area = area;
            k = 1;
            
            %Update bc
            ldir = bc.Neu.ldir;
            ldof = bc.Neu.ldof;
            lval = bc.Neu.lval;
            pbc = bouncon(p,e,t,ldir,ldof,lval);
            pbc = unique(pbc,'rows'); %just for symetry_cog_2 eliminate otherwise
            bc.pNeu = pbc;
            
            ldir = bc.Dir.ldir;
            ldof = bc.Dir.ldof;
            lval = bc.Dir.lval;
            pbc = bouncon(p,e,t,ldir,ldof,lval);
            pbc = unique(pbc,'rows'); %just for symetry_cog_2 eliminate otherwise
            bc.pDir = pbc;
            
            % freeze groups
%             phold = [];
%             for i = 1:size(mesh.ghold,1)
%                 phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i)))); %nodes to block
%             end
%             pfree = setdiff(1:size(p,2),phold);

phold = [];

rg = mesh.RG(:)';
ghold = mesh.ghold(:)';

for i = 1:length(ghold)
    phold = [phold; unique(t(1:3,rg == ghold(i)))];
end

phold = unique(phold);

pfree = setdiff(1:size(p,2),phold);

            option = 'null';
            
            % shape functional and volume associated to the hold-all domain
            psi_hold_all = ones(length(p),length(matprop.E)-1); psi_hold_all(:,1) = -1;
            [U,F,vol_hold_all] = solve(psi_hold_all, mesh, matprop, pdecoef, bc, signatures);
            params.energy0 = 0.5 * dot(F,U); % initial energy
            max_vol = vol_hold_all(1);
            params.volinit = max_vol; % store maximum volume for each fase
            voltarget = max_vol.*volfrac;
            
            % stop criterion over the volume constraint
            volstop = voleps.*max_vol;
            if penalization == 1
                volstop = max_vol;
            end
        end
    end    
end

fig = input(' -> export figures? "y" or "n"        : ', 's');
if strcmp(fig,'y')
    cd('results')
    save(strcat(filename,'.mat'))
        print (1, '-dtiff',  strcat(filename,'_meshfem'));
        print (2, '-dtiff', strcat(filename,'_iniguess'));
%         print (3, '-djpeg', strcat(filename,'_topology'));
        saveas(3,strcat(filename,'_topology','.fig'))
        print (4, '-dtiff',  strcat(filename,'_shapefunc'));
        print (5, '-dtiff',  strcat(filename,'_thetaangle'));
        print (6, '-dtiff',  strcat(filename,'_volume'));
    cd ..
else
    return;
end