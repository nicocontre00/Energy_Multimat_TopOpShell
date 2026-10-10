function [wing] = transformation

% Point sets 
p1 = [2.686041355;0.3291043639;3.116122484];
p2 = [2.685811281;0.3303097188;3.116122484]; 
p3 = [1.318630576;0.04338237271;0.3116122484];
r1 = [2.686041761; 0.2383747116; 3.116122402]; 
r2 = [2.685810076; 0.2395798007; 3.116122402]; 
r3 = [1.318637301; -0.006605711034; 0.3116122402];

p4 = [0.2449964434; 0.02215306275; 0.3116122484];
r4 = [0.2449964363; 0.02215306169; 0.3116122402];

p5 = [2.449964285; 0.2843421102; 3.116122484];
r5 = [2.449964363; 0.2843421111; 3.116122402];

p6 = [1.676238298; 0.1724174023; 1.253752351];
r6 = [1.679855542; 0.01009001659; 1.262666248];

p = [p1,p2,p3,p4,p5,p6];
r = [r1,r2,r3,r4,r5,r6];

% centroid of the sets
c_p = mean(p,2);
c_r = mean(r,2);

Q_p = [p1,p2,p3,p4,p5,p6]-c_p;
Q_r = [r1,r2,r3,r4,r5,r6]-c_r;

% H = Q_p * Q_r';
% 
% % Find the SVD of H (H = U * A * V')
% [U, A, V] = svd(H);
% X = V*U';
% 
% if det(X) > 0
% 	R = X; % Rotation matrix
% 	T = c_r - R*c_p; % translation matrix
% else
% 	R = [];
% end

H = Q_p * Q_r';
[U,~,V] = svd(H);

R = V*U';
T = c_r - R*c_p;

error_mm = 1000*vecnorm(R*p + T - r,2,1)

cd('examples\shell')
run('wing_barycenters.m');
g = wing.g;

% wing.g0 = (R*g')' + T'; % Transformed barycenters matrix

wing.g0 = g

end
