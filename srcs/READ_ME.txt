Running 'calibrate.m' it will be displayed on console the calibration parameter
taken from Xsens™ MTi™ we use on our real tests followed by the estimated parameters
we obtain with our method.

'IMU0x2Dalpha.mat' and 'IMU0x2Domega.mat' are the accelerometer raw data matrix and
the gyroscope raw data matrix respectively.

All the other file are functions used into 'calibrate.m':
- 'accCostFunctLSQNONLIN.m' is the cost function used to optimize accelerometer
	sensor error model parameters;
- 'gyroCostFunctLSQNONLIN.m' is the cost function used to optimize gyroscope
	sensor error model parameters;
- 'rotationRK4.m' is the implementation of the quaternion kinematics Runge-Kutta
	4th order integration algorithm;
- 'fromOmegaToQ. m' compute the unit quaternion q associated to a three-dimensional 
	angular velocity and a period of time t;
- 'fromQtoR.m' compute the rotation matrix associated to a unit quaternion q;
- 'obtainComparabelMatrix.m' compute the transformation that permit to compare the
	estimated matices to the ones given in the datasheet.