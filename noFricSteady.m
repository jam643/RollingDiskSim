function noFricSteady
%use symbolic algebra to calculate circular motion for no friction
syms phi theta psi wp ws 'real'
syms R m g 'real'

i=[1,0,0]'; j=[0,1,0]'; k=[0,0,1]';

et=-sin(phi)*i+cos(phi)*j; 
ip=cos(phi)*i+sin(phi)*j;
lambda=cos(theta)*k+sin(theta)*ip; 
n=cos(theta)*ip-sin(theta)*k;

etd=cross(wp*k,et);
lambdad=cross(wp*k,lambda);
nd=cross(wp*k,n);

rGrelC=R*lambda;

w=wp*k+ws*n;
alpha=ws*nd;

vG=cross(w,R*lambda);
aG=0*i+0*j+0*k;

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

x=solve(dot(AMB,i),wp);
eq=simplify(x)