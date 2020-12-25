function zdot=noSlipDisk_rhs(t,z,p)

m=p.m; R=p.R; g=p.g;                                  
phi    = z(1)             ;                          
theta = z(2)             ;                          
psi    = z(3)             ;                          
phid = z(4)             ;                        
thetad = z(5)             ;                        
psid = z(6)             ;                        
xG = z(7)             ;                        
yG = z(8)             ;                        

phidd = ...
(2*psid*thetad)/cos(theta);                         

thetadd = ...
(4*g*sin(theta) + 5*R*phid^2*cos(theta)*sin(theta) - 6*R*phid*psid*cos(theta))/(5*R);                         

psidd = ...
(thetad*(5*phid*cos(theta)^2 + 6*psid*sin(theta)))/(3*cos(theta));                         

xGdot = ...
R*cos(theta)*(thetad*cos(phi) + psid*cos(theta)*sin(phi)) - R*sin(phi)*sin(theta)*(phid - psid*sin(theta));                         

yGdot = ...
R*cos(theta)*(thetad*sin(phi) - psid*cos(phi)*cos(theta)) + R*cos(phi)*sin(theta)*(phid - psid*sin(theta));                         

zdot = [phid thetad psid phidd thetadd psidd xGdot yGdot]'  ;            
