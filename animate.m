function animate(tarray,zarray,p)
%%%animates no slip and no friction disks

R=p.R;
N=30;
disk=zeros(3,N);
%create disk
diskRot=disk;
for k=1:N+1
    disk(2:3,k)=R*[cos(2*pi*(k-1)/N),sin(2*pi*(k-1)/N)]';
end
%create spokes
spoke1=R*[0,0,0,0;1,cos(2*pi/3),cos(4*pi/3),1;0,sin(2*pi/3),sin(4*pi/3),0];
spoke2=0.5*R*[0,0,0,0;cos(pi/3),-1,cos(5*pi/3),cos(pi/3);sin(pi/3),0,sin(5*pi/3),sin(pi/3)];

%create axis limits
axisLim=[min(zarray(:,7))-R,max(zarray(:,7))+R,min(zarray(:,8))-R,max(zarray(:,8))+R,0,2*R];

%animation
f=figure;
set(f,'units','normalized','outerposition',[0 0 1 1],'color','w');
for t=1:length(tarray)
    if ~ishandle(f)
        break;
    end
    %unpack/calculate variables
    phi=zarray(t,1); theta=zarray(t,2); psi=zarray(t,3);
    xG=zarray(t,7); yG=zarray(t,8);
    %calculate rotation matrix
    Rotz=[cos(phi), -sin(phi),  0;...
        sin(phi),  cos(phi),  0;...
        0,          0,      1];
    Roty=[cos(theta),   0,  sin(theta);...
            0,          1,      0;...
        -sin(theta),    0,   cos(theta)];
    Rotx=[1,   0,  0;...
            0, cos(psi),      -sin(psi);...
        0,    sin(psi),   cos(psi)];
    Rot=Rotz*Roty*Rotx;
    %rotate disk
    for k=1:size(disk,2)
        diskRot(:,k)=Rot*disk(:,k)+[xG,yG,R*cos(theta)]';
    end
    %rotate spokes
    for k=1:size(spoke1,2)
        spoke1Rot(:,k)=Rot*spoke1(:,k)+[xG,yG,R*cos(theta)]';
        spoke2Rot(:,k)=Rot*spoke2(:,k)+[xG,yG,R*cos(theta)]';
    end
    %plot
    cla;
    fill3(diskRot(1,:),diskRot(2,:),diskRot(3,:),[0.9,0.9,0.9],'linewidth',3);
    hold on;
    plot3(spoke1Rot(1,:),spoke1Rot(2,:),spoke1Rot(3,:),'k','linewidth',3);
    plot3(spoke2Rot(1,:),spoke2Rot(2,:),spoke2Rot(3,:),'k','linewidth',3);
    
    %plot point C over time
    ip=cos(phi)*[1,0,0]'+sin(phi)*[0,1,0]';
    lambda=cos(theta)*[0,0,1]'+sin(theta)*ip;
    xC(t)=xG-R*lambda(1);
    yC(t)=yG-R*lambda(2);    
    plot(xC(1:t),yC(1:t),'k','linewidth',3);
    plot3(zarray(1:t,7),zarray(1:t,8),R*cos(zarray(1:t,2)),'r','linewidth',1);
    
    axis equal;
    axis(axisLim);
    pause(0.001);
end