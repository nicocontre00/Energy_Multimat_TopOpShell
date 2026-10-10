function [rg_p] = rgpernode(p,t,r,g)

num_p = size(p,2);
rg_p = []; 

for j = 1:num_p	
	[~, idj] = find(t(:,1:3) == j);     
	rg_j = unique(r(:, idj));	
	if length(rg_j) > 1		
		p_rgj = p(j,:);
		g_rgj = g(ismember(g, rg_j),:);		
		[~,id] = min((p_rgj(1) - g_rgj(1,:)).^2 + (p_rgj(2) - g_rgj(2,:)).^2 + (p_rgj(3) - g_rgj(3,:)).^2);		       
		rgj = g_rgj(id,:);	
		rgp = [rg_p ; rgj];
    else
		rgp = [rg_p ; rg_j]; 
    end	
end