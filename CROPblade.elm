# CROP simple MBDyn model, Louis Gagnon 2014
set: CURR_BLADE_ANGLE = 2. * pi * (CURR_BLADE/1000.-1.) / NBLADES ; # angle from positive x-axis, with positive dir. counter-clockwise
set: ti = -CURR_BLADE_ANGLE / OMEGA;


joint: CURR_DRUM + CURR_BLADE + 1, revolute hinge,  # blade free rotation
	 CURR_DRUM,
	 position, reference, node, ROTOR_DIAMETER/2. * cos(CURR_BLADE_ANGLE), ROTOR_DIAMETER/2. * sin(CURR_BLADE_ANGLE), 0,
	orientation, 
		1, 1., 0., 0.,
		3, 0., 0., 1.,
	 CURR_DRUM + CURR_BLADE,
	 position, reference, node, 0,0,0,
	orientation, 
		1, 1., 0., 0.,
		3, 0., 0., 1.;
	
joint: CURR_DRUM + CURR_BLADE + 2, total joint,  # blade imposed rotation
	 CURR_DRUM,
	 position, reference, node, ROTOR_DIAMETER/2. * cos(CURR_BLADE_ANGLE), ROTOR_DIAMETER/2. * sin(CURR_BLADE_ANGLE), 0,
	 CURR_DRUM + CURR_BLADE,
	 position, reference, node, 0,0,0,
	position constraint,
		0, 0, 0,
		null,
	orientation constraint,
		0, 0, 1,
		0, 0, 1,
#		CURR_BLADE_ANGLE;
	array,
	5.,
	theta_o + CURR_BLADE_ANGLE,  # constant drive (collective pitch + initial tangential position of blade to drum)
	cosine, ti , OMEGA, -theta_c, forever, theta_c,
	sine, ti , OMEGA, theta_s, forever, 0.,
	cosine, ti , 2.*OMEGA, theta_c2, forever, 0., ### TODO: fix as above
	sine, ti , 2.*OMEGA, theta_s2, forever, 0.;
	
	remark: "ELM: blade, CURR_BLADE_ANGLE, INIT_OMG_OSC", CURR_BLADE/BLADE_1, CURR_BLADE_ANGLE/pi*180, INIT_OMG_OSC / pi*180 ;
	remark: "cos(pi)", cos(pi);
	
aerodynamic body: CURR_DRUM + CURR_BLADE,
	CURR_DRUM + CURR_BLADE,
	user defined induced velocity, CURR_INFLOW,
	reference, node,
	0., 0., 0.,
	reference, node,
		1, 0., 1., 0.,
		2, 1., 0., 0.,
	SPAN,
	const, CHORD,	# chord
	const, CHORD/4,	# aero center
	const, BC_OFFSET,	# b c point
	const, 0.,	# twist
#	 tip loss, const, 0.95,
	AERO_INTEGRATION_POINTS,
	# theodorsen,
		c81, NACA0012,
	jacobian, no;
#	jacobian, AERO_JACOBIAN;
	
	
	
