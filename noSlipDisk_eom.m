function noSlipDisk_eom()
%calculates symbolically equations of motion for a rolling disk

syms phi theta psi phid thetad psid phidd thetadd psidd 'real'
syms R m g 'real'

%unit vectors
i=[1,0,0]'; j=[0,1,0]'; k=[0,0,1]';

et=-sin(phi)*i+cos(phi)*j; 
ip=cos(phi)*i+sin(phi)*j;
lambda=cos(theta)*k+sin(theta)*ip; 
n=cos(theta)*ip-sin(theta)*k;

%rate of change or unit vectors
etd=cross(phid*k,et);
lambdad=cross((phid*k+thetad*et),lambda);
nd=cross(phid*k+thetad*et,n);

%position vector
rGrelC=R*lambda;

%angular velocity and acceleration
w=phid*k+thetad*et+psid*n;
alpha=phidd*k+thetadd*et+thetad*etd+psidd*n+psid*nd;

%velocity and acceleration
vG=cross(w,R*lambda);
aG=cross(alpha,R*lambda)+cross(w,R*lambdad);

%rotation matrix
Rotz=[cos(phi), -sin(phi),  0;...
    sin(phi),  cos(phi),  0;...
    0,          0,      1];
Roty=[cos(theta),   0,  sin(theta);...
    0,          1,      0;...
    -sin(theta),    0,   cos(theta)];

Rot=Rotz*Roty;

%moment of inertia of uniform disk converted to fixed frame
I=m*R^2/4;
IGrelB=[2*I,0,0;0,I,0;0,0,I];
IGrelF=Rot*IGrelB*Rot.';

%AMB
HdotrelC=cross(rGrelC,m*aG)+IGrelF*alpha+cross(w,IGrelF*w);

MrelC=cross(rGrelC,-m*g*k);

AMB=MrelC-HdotrelC;

%solve for the double dots
A=jacobian(AMB,[phidd,thetadd,psidd]');
B=AMB-A*[phidd,thetadd,psidd]';
x=simplify(-(A^-1)*B);

%write rhs file containing phidd, thetadd, psidd for ode45
fid=fopen(   'noSlipDisk_rhs.m','w'                                    );
fprintf(fid, 'function zdot=noSlipDisk_rhs(t,z,p)\n\n');
fprintf(fid, 'm=p.m; R=p.R; g=p.g;                                  \n');
fprintf(fid, 'phi    = z(1)             ;                          \n');
fprintf(fid, 'theta = z(2)             ;                          \n');
fprintf(fid, 'psi    = z(3)             ;                          \n');
fprintf(fid, 'phid = z(4)             ;                        \n');
fprintf(fid, 'thetad = z(5)             ;                        \n');
fprintf(fid, 'psid = z(6)             ;                        \n');
fprintf(fid, 'xG = z(7)             ;                        \n');
fprintf(fid, 'yG = z(8)             ;                        \n\n');

fprintf(fid,['phidd = ...\n' char(x(1)) ';                         \n\n']);
fprintf(fid,['thetadd = ...\n' char(x(2)) ';                         \n\n']);
fprintf(fid,['psidd = ...\n' char(x(3)) ';                         \n\n']);
fprintf(fid,['xGdot = ...\n' char(vG(1)) ';                         \n\n']);
fprintf(fid,['yGdot = ...\n' char(vG(2)) ';                         \n\n']);

fprintf(fid, 'zdot = [phid thetad psid phidd thetadd psidd xGdot yGdot]''  ;            \n');
fclose(fid);