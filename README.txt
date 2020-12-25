6700 no slip/ no friction disk in 3D final project
Jesse Miller jam643

To simulate rolling disk with no slip:
-type ROOTnoSlip.m in command window to run, animate, and plot
-set circle=0 to be able to input initial conditions
	-can change these initial conditions:
	phi theta psi phid thetad psid xG yG
-set circle=1 for circular motion conditions
	-can change these initial conditions:
	theta(between -pi/2 to pi/2) ws(as long as it's greater than the minimum spin rate wsmin)
-set circle=2 for intersection with no friction circular motion
	-can only change theta (between -pi/2 to pi/2)

To simulate slidding disk with no friction:
-type ROOTnoFric.m in command window to run, animate, and plot
-set circle=0 to be able to input initial conditions
	-can change these initial conditions:
	phi theta psi phid thetad psid xG yG xGd yGd
-set circle=1 for circular motion conditions
	-can change these initial conditions:
	theta(between -pi/2 to pi/2) ws(as long as it's greater than the minimum spin rate wsmin)
-set circle=2 for intersection with no slip circular motion
	-can only change theta (between -pi/2 to pi/2)