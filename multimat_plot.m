function multimat_plot( p,e,t,fi,AX )
%multimat_plot Plots the caracteristic function fi of each material and
%also the void region. Plots up to 4 different materials and void.

%auth: Augusto Romero
%date: 26-06-2018
%updated: 09-06-2025 - Matlab does not suppor anymore pdesurf function

n = size(fi,2); np = max(size(p));
colors = [0 0 0;1 0 0; 0 0 1;0 1 0];
PATCHES = findall(AX,'type','patch');
if or(isempty(PATCHES), np~=n) 
    axes(AX);
    for i=1:n	
		pdeplot(p,e,t,'XYdata',double(fi(:,i)),'ZData',double(fi(:,i)),'ColorBar','off');
 %      pdeplot(p,e,t,'XYdata',double(fi(:,i)),'ZData',double(fi(:,i)),'ColorBar','off','Parent',AX);
        set(AX, 'NextPlot', 'add');
        PATCHES = findall(AX,'type','patch');
        PATCHES(1).FaceColor = colors(i,:);
    end
    AX.View = [0 90]; AX.DataAspectRatio = [1 1 1];
    AX.ZLim = [1e-1 1.2];
    axis off; 
    
    PATCHES(1).FaceColor = [1 1 1];
else
    for i=1:n
        PATCHES(n+1-i).Vertices(:,3)=fi(:,i);
    end
end



