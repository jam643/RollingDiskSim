%uses symbolic notation to calculate the parameters for circular motion for
%rolling disk
syms phi theta psi ws 'real'
syms R m g I 'real'

wp=ws/sin(theta);

i=[1,0,0]'; j=[0,1,0]'; k=[0,0,1]';

et=-sin(phi)*i+cos(phi)*j; 
ip=cos(phi)*i+sin(phi)*j;
lambda=cos(theta)*k+sin(theta)*ip; 
n=cos(theta)*ip-sin(theta)*k;

etd=cross(wp*k,et);
lambdad=cross((wp*k),lambda);
nd=cross(wp*k,n);

rGrelC=R*lambda;

w=wp*k+ws*n;
alpha=ws*nd;

vG=cross(w,R*lambda);
aG=cross(alpha,R*lambda)+cross(w,R*lambdad);

Rotz=[cos(phi), -sin(phi),  0;...
    sin(phi),  cos(phi),  0;...
    0,          0,      1];
Roty=[cos(theta),   0,  sin(theta);...
    0,          1,      0;...
    -sin(theta),    0,   cos(theta)];

Rot=Rotz*Roty;

I=m*R^2/4;
IGrelB=[2*I,0,0;0,I,0;0,0,I];
IGrelF=Rot*IGrelB*Rot.';

HdotrelC=cross(rGrelC,m*aG)+IGrelF*alpha+cross(w,IGrelF*w);

MrelC=cross(rGrelC,-m*g*k);

AMB=MrelC-HdotrelC;

x=solve(dot(AMB,i),ws);
eq=simplify(x)


