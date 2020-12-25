function plotNoSlip(tarray,zarray,p)
%create plots for no slip disk

f=figure;
set(f,'units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1);
plot(tarray,zarray(:,1),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('Phi [rad]','fontsize',14);
subplot(2,3,2);
plot(tarray,zarray(:,2),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('Theta [rad]','fontsize',14);
subplot(2,3,3);
plot(tarray,zarray(:,3),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('Psi [rad]','fontsize',14);
subplot(2,3,4);
plot(tarray,zarray(:,4),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('Phi dot [rad/s]','fontsize',14);
subplot(2,3,5);
plot(tarray,zarray(:,5),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('Theta dot [rad/s]','fontsize',14);
subplot(2,3,6);
plot(tarray,zarray(:,6),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('Psi dot [rad/s]','fontsize',14);

f=figure;
set(f,'units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(tarray,zarray(:,7),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('x [m]','fontsize',14);
subplot(2,2,2);
plot(tarray,zarray(:,8),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('y [m]','fontsize',14);

for t=1:length(tarray)
    phi=zarray(t,1); theta=zarray(t,2); psi=zarray(t,3);
    phid=zarray(t,4); thetad=zarray(t,5); psid=zarray(t,6);
    m=p.m; R=p.R;
    xG=zarray(t,7); yG=zarray(t,8);
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
    i=[1,0,0]'; j=[0,1,0]'; k=[0,0,1]';
    
    et=-sin(phi)*i+cos(phi)*j;
    ip=cos(phi)*i+sin(phi)*j;
    lambda=cos(theta)*k+sin(theta)*ip;
    n=cos(theta)*ip-sin(theta)*k;
    
    w=phid*k+thetad*et+psid*n;
    
    vG(t,1:3)=cross(w,R*lambda);
    
    E(t)=0.5*m*dot(vG(t,:),vG(t,:))+0.5*dot(w,(IGrelF*w))+m*p.g*R*cos(theta);
end

subplot(2,2,3);
plot(tarray,vG(:,1),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('x dot [m/s]','fontsize',14);
subplot(2,2,4);
plot(tarray,vG(:,2),'linewidth',2);
xlabel('Time [s]','fontsize',14);
ylabel('y dot [m/s]','fontsize',14);

f=figure;
set(f,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
plot(tarray,E);
xlabel('Time [s]','fontsize',14);
ylabel('Energy [Joule]','fontsize',14);
subplot(1,2,2)
plot(tarray,E-E(1));
xlabel('Time [s]','fontsize',14);
ylabel('Change in Energy [Joule]','fontsize',14);