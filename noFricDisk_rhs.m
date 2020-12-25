function zdot=noFricDisk_rhs(t,z,p)

m=p.m; R=p.R; g=p.g;                                  
phi    = z(1)             ;                          
theta = z(2)             ;                          
psi    = z(3)             ;                          
phid = z(4)             ;                        
thetad = z(5)             ;                        
psid = z(6)             ;                        
xG = z(7)             ;                        
yG = z(8)             ;                        
xGd = z(9)             ;                        
yGd = z(10)             ;                        

phidd = ...
(2*psid*thetad)/cos(theta);                         

thetadd = ...
-(8*g*sin(theta) + R*phid^2*sin(2*theta) - 4*R*thetad^2*sin(2*theta) - 4*R*phid*psid*cos(theta))/(2*R*(4*cos(theta)^2 - 5));                         

psidd = ...
(thetad*(phid*cos(theta)^2 + 2*psid*sin(theta)))/cos(theta);                         

xGdd =0;                        

yGdd =0;                     

zdot = [phid thetad psid phidd thetadd psidd xGd yGd xGdd yGdd]'  ;            
