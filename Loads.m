% clear all, close all, clc
function [pNeu] = Loads

cd('examples\shell')
run('wing_msh.m');
p = msh.POS; t = msh.TRIANGLES; r = msh.RG;

num_r = msh.shell; % Geometric region number of the new mesh
% t_r = t(ismember(r,num_r),1:3);
% id1 = t_r(:,1); id2 = t_r(:,2); id3 = t_r(:,3);
% x = [p(id1,1), p(id2,1), p(id3,1)];
% y = [p(id1,2), p(id2,2), p(id3,2)];
% z = [p(id1,3), p(id2,3), p(id3,3)];
% g = [mean(x,2), mean(y,2), mean(z,2)]; % New barycenters matrix

top = msh.topface; bot = msh.bottomface; side = msh.sideface;
ttop = t(ismember(r,top),1:3); tbot = t(ismember(r,bot),1:3); tside = t(ismember(r,side),1:3);
idtop1 = ttop(:,1); idtop2 = ttop(:,2); idtop3 = ttop(:,3); idtop = [idtop1;idtop2;idtop3];
idbot1 = tbot(:,1); idbot2 = tbot(:,2); idbot3 = tbot(:,3); idbot = [idbot1;idbot2;idbot3];
idside1 = tside(:,1); idside2 = tside(:,2); idside3 = tside(:,3); idside = [idside1;idside2;idside3];
gtop = (p(ttop(:,1),:)  + p(ttop(:,2),:)  + p(ttop(:,3),:))/3; % Barycenters matrix of the top face
gbot = (p(tbot(:,1),:)  + p(tbot(:,2),:)  + p(tbot(:,3),:))/3; % Barycenters matrix of the bottom face
gside = (p(tside(:,1),:) + p(tside(:,2),:) + p(tside(:,3),:))/3; % Barycenters matrix of the side face

run('wing_oldmsh.m');
p_old = msh.POS_old; t_old = msh.TRIANGLES_old; r_old = msh.RG_old;
top0 = msh.topface0; bot0 = msh.bottomface0; side0 = msh.sideface0;

id1_old = t_old(:,1); id2_old = t_old(:,2); id3_old = t_old(:,3);
x_old = [p_old(id1_old,1), p_old(id2_old,1), p_old(id3_old,1)];
y_old = [p_old(id1_old,2), p_old(id2_old,2), p_old(id3_old,2)];
z_old = [p_old(id1_old,3), p_old(id2_old,3), p_old(id3_old,3)];
g_old = [mean(x_old,2), mean(y_old,2), mean(z_old,2)]; % Old barycenters matrix
g0top = g_old(ismember(r_old,top0),:); % Old barycenters matrix of the top face
g0bot = g_old(ismember(r_old,bot0),:); % Old barycenters matrix of the bottom face
g0side = g_old(ismember(r_old,side0),:); % Old barycenters matrix of the side face

run('wing_coefficients.m'); 
Cf0_t = wing.Cft / 150; Cf0_z = wing.Cfz / 300; Cp0 = wing.Cp;  % Old Wing Coefficients
% Cf0_t = pdeprtni(p_old, t_old, Cf0_t);
% Cf0_z = pdeprtni(p_old, t_old, Cf0_z);
% Cp0 = pdeprtni(p_old, t_old, Cp0);

Cp0top = Cp0(ismember(r_old,top0),:); % Old wing Coefficients of the top face
Cp0bot = Cp0(ismember(r_old,bot0),:); % Old wing Coefficients of the bottom face
Cp0side = Cp0(ismember(r_old,side0),:); % Old wing Coefficients of the side face
Cf0_ztop = Cf0_z(ismember(r_old,top0),:); % Old wing Coefficients of the top face
Cf0_zbot = Cf0_z(ismember(r_old,bot0),:); % Old wing Coefficients of the bottom face
Cf0_zside = Cf0_z(ismember(r_old,side0),:); % Old wing Coefficients of the side face
Cf0_ttop = Cf0_t(ismember(r_old,top0),:); % Old wing Coefficients of the top face
Cf0_tbot = Cf0_t(ismember(r_old,bot0),:); % Old wing Coefficients of the bottom face
Cf0_tside = Cf0_t(ismember(r_old,side0),:); % Old wing Coefficients of the side face

ICptop = scatteredInterpolant(g0top(:,1),g0top(:,2),g0top(:,3),Cp0top, 'natural', 'nearest');
ICpbot = scatteredInterpolant(g0bot(:,1),g0bot(:,2),g0bot(:,3),Cp0bot, 'natural', 'nearest');
ICpside = scatteredInterpolant(g0side(:,1),g0side(:,2),g0side(:,3),Cp0side, 'natural', 'nearest');
ICf_ztop = scatteredInterpolant(g0top(:,1),g0top(:,2),g0top(:,3),Cf0_ztop, 'natural', 'nearest');
ICf_zbot = scatteredInterpolant(g0bot(:,1),g0bot(:,2),g0bot(:,3),Cf0_zbot, 'natural', 'nearest');
ICf_zside = scatteredInterpolant(g0side(:,1),g0side(:,2),g0side(:,3),Cf0_zside, 'natural', 'nearest');
ICf_ttop = scatteredInterpolant(g0top(:,1),g0top(:,2),g0top(:,3),Cf0_ttop, 'natural', 'nearest');
ICf_tbot = scatteredInterpolant(g0bot(:,1),g0bot(:,2),g0bot(:,3),Cf0_tbot, 'natural', 'nearest');
ICf_tside = scatteredInterpolant(g0side(:,1),g0side(:,2),g0side(:,3),Cf0_tside, 'natural', 'nearest');

Cptop = ICptop(gtop(:,1),gtop(:,2),gtop(:,3));
Cpbot = ICpbot(gbot(:,1),gbot(:,2),gbot(:,3));
Cpside = ICpside(gside(:,1),gside(:,2),gside(:,3));
Cf_ztop = ICf_ztop(gtop(:,1),gtop(:,2),gtop(:,3));
Cf_zbot = ICf_zbot(gbot(:,1),gbot(:,2),gbot(:,3));
Cf_zside = ICf_zside(gside(:,1),gside(:,2),gside(:,3));
Cf_ttop = ICf_ttop(gtop(:,1),gtop(:,2),gtop(:,3));
Cf_tbot = ICf_tbot(gbot(:,1),gbot(:,2),gbot(:,3));
Cf_tside = ICf_tside(gside(:,1),gside(:,2),gside(:,3));

Cp = [Cptop;Cpbot;Cpside];
Cf_z = [Cf_ztop;Cf_zbot;Cf_zside];
Cf_t = [Cf_ttop;Cf_tbot;Cf_tside];

% cd ..\..
% [rgp] = rgpernode(p,t,r,g); % Geometric region per node of the new mesh
% [rgp_old] = rgpernode(p_old,t_old,r_old,g_old); % Geometric region per node of the new mesh

% [wing] = tranformation; % Old barycenters matrix
% g0 = wing.g0; g0 = g0*1000;
% 
% ICp = scatteredInterpolant(g0(:,1),g0(:,2),g0(:,3),Cp0);
% ICf_z = scatteredInterpolant(g0(:,1),g0(:,2),g0(:,3),Cf0_z);
% ICf_t = scatteredInterpolant(g0(:,1),g0(:,2),g0(:,3),Cf0_t);
% Cp = ICp(g(:,1),g(:,2),g(:,3)); % New Wing Coefficients
% Cf_t = ICf_t(g(:,1),g(:,2),g(:,3)); % New Wing Coefficients
% Cf_z = ICf_z(g(:,1),g(:,2),g(:,3)); % New Wing Coefficients

t_r = [ttop;tbot;tside];
id1 = t_r(:,1);
id2 = t_r(:,2);
id3 = t_r(:,3);
v1 = p(id1,:); v2 = p(id2,:); v3 = p(id3,:); 
n = cross(v2-v1,v3-v1,2); 
A = vecnorm(n,2,2)/2;
n = n ./ vecnorm(n,2,2);

itop = (1:size(ttop,1))';
ibot = size(ttop,1) + (1:size(tbot,1))';
invertir = [itop(n(itop,2)<0); ibot(n(ibot,2)>0)];
n(invertir,:) = -n(invertir,:);

iside = size(ttop,1) + size(tbot,1) + (1:size(tside,1))';
invertir = iside(n(iside,1)<0);
n(invertir,:) = -n(invertir,:);

q = 0.0512084; % [MPa] 51208.4 Pa
Fp = -q.*A.*Cp.*n;

Cf = [Cf_t, zeros(size(Cf_t,1),1), Cf_z];
Ff = q.*A.*(Cf - dot(Cf,n,2).*n);
F = Fp + Ff;

node = t_r'; node = node(:); node = [node; node; node];
f = repelem(F/3,3,1); f = [f(:,1); f(:,2); f(:,3)]; % [fx,fy,fz]
dof1 = ones(size(f,1)/3,1); dof2 = 2*ones(size(f,1)/3,1); dof3 = 3*ones(size(f,1)/3,1); 
dof = [dof1;dof2;dof3];
Neu = [node, dof, f];

[dofs,~,idx] = unique(Neu(:,1:2),'rows');
pNeu = [dofs, accumarray(idx,Neu(:,3))];

% Neu = [];
% for i = 1 : size(pNeu,1)   
%     row = pNeu(i,:);    
%     idx = pNeu(:,1) == row(1,1) & pNeu(:,2) == row(1,2);
% 
%     if isempty(Neu)
%         Neu(1,:) = [row(1,1), row(1,2), sum(pNeu(idx,3))];
%     else
%         if sum(Neu(:,1) == row(1,1) & Neu(:,2) == row(1,2))==0
%             Neu(end+1,:) = [row(1,1), row(1,2), sum(pNeu(idx,3))]; 
%         end
%     end
% end
% 
% pNeu = Neu;