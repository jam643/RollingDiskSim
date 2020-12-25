function ROOTnoSlip()
%calculates and animates a disk sliding on a surface with no friction
close all;
%parameters
p.m=1; p.g=10; p.R=1;

circle=0; %0: select any parameters 1: only circular motion 2: intersection with no slip
if circle==2
    %intersection between no friction circular motion
    g=p.g; R=p.R;
    theta=pi/3; %can only edit this value (-pi/2:pi/2) for precession motion
    phi0=0; psi0=0; thetad0=0;
    %spin rate for intersection with no friction
    ws=(2*sin(theta)*(R*g*cos(theta))^(1/2))/(R*cos(theta));
    %rate of precession
    wp=(((R*cos(theta)*(20*g*cos(theta)^2 - 20*g + 9*R*ws^2*cos(theta)))^(1/2) + 3*R*ws*cos(theta))/(5*R*cos(theta)*sin(theta)));
    phid0=wp; psid0=ws; theta0=theta;
    %analytic radius
    radAnalytic=R*abs((sin(theta)-ws/wp));
elseif circle==1
    %circular motion initial conditions
    g=p.g; R=p.R;
    theta=pi/3; %can only edit this value (-pi/2:pi/2) for precession motion
    phi0=0; psi0=0; thetad0=0;
    %minimum spin rate for circular motion
    wsmin=sqrt(20*g/(9*R)*(1-cos(theta)^2)/cos(theta));
    %can set multiple of minimum spin rate
    ws=wsmin;
    %precession rate for circular motion
    wp=(((R*cos(theta)*(20*g*cos(theta)^2 - 20*g + 9*R*ws^2*cos(theta)))^(1/2) + 3*R*ws*cos(theta))/(5*R*cos(theta)*sin(theta)));
    phid0=wp; psid0=ws; theta0=theta;
    %analytic radius
    radAnalytic=R*abs((sin(theta)-ws/wp));    
else
    %initial conditions
    phi0=0; theta0=0; psi0=0;
    phid0=1; thetad0=2; psid0=3;
end
xG0=0; yG0=0;

%time vector
t=10; fps=50;
tspan=linspace(0,t,t*fps);

%store initial conditions
z0=[phi0, theta0, psi0, phid0, thetad0, psid0, xG0, yG0]';

options=odeset('abstol',1e-8,'reltol',1e-8);

%calculate and solve equations of motion
% noSlipDisk_eom();
[tarray,zarray]=ode45(@noSlipDisk_rhs,tspan,z0,options,p);

%displays difference in radius
if circle
    radNumeric=(max(zarray(:,7))-min(zarray(:,7)))/2;
    fprintf('Difference between radius and predicted = %0.2e meters\n',abs(radNumeric-radAnalytic));
end

%plot and animate
animate(tarray,zarray,p)
plotNoSlip(tarray,zarray,p);

