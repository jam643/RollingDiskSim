function ROOTnoFric()
%calculates and animates a disk sliding on a surface with no friction

close all;
%parameters
p.m=1; p.g=10; p.R=1;

circle=0; %0: select any parameters 1: only circular motion 2: intersection with no slip
if circle==2
    %chose conditions to create circular motion that overlap with no slip
    g=p.g; R=p.R;
    theta=pi/3;  %can only edit this value (-pi/2:pi/2) for precession motion
    phi0=0; psi0=0; thetad0=0;
    %spin rate for intersection with no slip
    ws=(2*sin(theta)*(R*g*cos(theta))^(1/2))/(R*cos(theta));
    %precession rate for circular motion
    wp=((R*cos(theta)*(4*g*cos(theta)^2 - 4*g + R*ws^2*cos(theta)))^(1/2) + R*ws*cos(theta))/(R*cos(theta)*sin(theta));
    phid0=wp; psid0=ws; theta0=theta;
    %analytic precession radius
    radAnalytic=0;
    xGd0=0; yGd0=0;
elseif circle==1
    %calculates circular motion conditions
    g=p.g; R=p.R;
    theta=pi/3; %can only edit this value (-pi/2:pi/2) for precession motion
    phi0=0; psi0=0; thetad0=0;
    %minimum spin for circular motion
    wsmin=sqrt(4*g*(1-cos(theta)^2)/(R*cos(theta)));
    %spin rate can be any value greater than wsmin
    ws=wsmin+0.0001;
    %precession rate for circular motion
    wp=((R*cos(theta)*(4*g*cos(theta)^2 - 4*g + R*ws^2*cos(theta)))^(1/2) + R*ws*cos(theta))/(R*cos(theta)*sin(theta));
    phid0=wp; psid0=ws; theta0=theta;
    radAnalytic=0;
    xGd0=0; yGd0=0;
else
    %initial conditions
    phi0=0; theta0=0; psi0=0;
    phid0=3; thetad0=1; psid0=10;
    xGd0=0; yGd0=0;
end
xG0=0; yG0=0;

%time vector
t=10; fps=50;
tspan=linspace(0,t,t*fps);

%store initial conditions
z0=[phi0, theta0, psi0, phid0, thetad0, psid0, xG0, yG0, xGd0, yGd0]';

options=odeset('abstol',1e-8,'reltol',1e-8);

%calculate and solve equations of motion
% noFricDisk_eom();
[tarray,zarray]=ode45(@noFricDisk_rhs,tspan,z0,options,p);

%display difference in calculated vs actual radius of precession
if circle
    radNumeric=(max(zarray(:,7))-min(zarray(:,7)))/2;
    fprintf('Difference between radius and predicted = %0.2e meters\n',abs(radNumeric-radAnalytic));
end

%animate and plot
animate(tarray,zarray,p)
plotNoFric(tarray,zarray,p)

