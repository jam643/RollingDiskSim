function noFricDisk_eom()
%calculates symbolically equations of motion for a rolling disk

syms phi theta psi phid thetad psid phidd thetadd psidd 'real'
syms R m g 'real'

i=[1,0,0]'; j=[0,1,0]'; k=[0,0,1]';

%unit vectors
et=-sin(phi)*i+cos(phi)*j; 
ip=cos(phi)*i+sin(phi)*j;
lambda=cos(theta)*k+sin(theta)*ip; 
n=cos(theta)*ip-sin(theta)*k;

%unit vector dots
etd=cross(phid*k,et);
lambdad=cross((phid*k+thetad*et),lambda);
nd=cross(phid*k+thetad*et,n);

%position vector
rGrelC=R*lambda;

%angular velocity and acceleration
w=phid*k+thetad*et+psid*n;
alpha=phidd*k+thetadd*et+thetad*etd+psidd*n+psid*nd;

%velocity and acceleration
aG=0*i+0*j+R*(-thetadd*sin(theta)-thetad^2*cos(theta))*k;

%Rotation matrix in banked frame
Rotz=[cos(phi), -sin(phi),  0;...
    sin(phi),  cos(phi),  0;...
    0,          0,      1];
Roty=[cos(theta),   0,  sin(theta);...
    0,          1,      0;...
    -sin(theta),    0,   cos(theta)];

Rot=Rotz*Roty;

%moment of inertia relative to body frame
I=m*R^2/4;
IGrelB=[2*I,0,0;0,I,0;0,0,I];
%moment of inertia relative to fixed frame
IGrelF=Rot*IGrelB*Rot.';

%AMB
HdotrelC=cross(rGrelC,m*aG)+IGrelF*alpha+cross(w,IGrelF*w);

MrelC=cross(rGrelC,-m*g*k);

AMB=MrelC-HdotrelC;

%solve for the double dots using jacobian
A=jacobian(AMB,[phidd,thetadd,psidd]');
B=AMB-A*[phidd,thetadd,psidd]';
x=simplify(-(A^-1)*B);

%write rhs file containing phidd, thetadd, psidd for ode45
fid=fopen(   'noFricDisk_rhs.m','w'                                    );
fprintf(fid, 'function zdot=noFricDisk_rhs(t,z,p)\n\n');
fprintf(fid, 'm=p.m; R=p.R; g=p.g;                                  \n');
fprintf(fid, 'phi    = z(1)             ;                          \n');
fprintf(fid, 'theta = z(2)             ;                          \n');
fprintf(fid, 'psi    = z(3)             ;                          \n');
fprintf(fid, 'phid = z(4)             ;                        \n');
fprintf(fid, 'thetad = z(5)             ;                        \n');
fprintf(fid, 'psid = z(6)             ;                        \n');
fprintf(fid, 'xG = z(7)             ;                        \n');
fprintf(fid, 'yG = z(8)             ;                        \n');
fprintf(fid, 'xGd = z(9)             ;                        \n');
fprintf(fid, 'yGd = z(10)             ;                        \n\n');

fprintf(fid,['phidd = ...\n' char(x(1)) ';                         \n\n']);
fprintf(fid,['thetadd = ...\n' char(x(2)) ';                         \n\n']);
fprintf(fid,['psidd = ...\n' char(x(3)) ';                         \n\n']);
fprintf(fid,'xGdd =0;                        \n\n');
fprintf(fid,'yGdd =0;                     \n\n');

fprintf(fid, 'zdot = [phid thetad psid phidd thetadd psidd xGd yGd xGdd yGdd]''  ;            \n');
fclose(fid);