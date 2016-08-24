from __future__ import division;

class usc:
	rc = 5.5;
	mO = 15.9994;
	mH = 1.00794;
	mSi = 28.0855;
	m = [mSi, mO, mH];
	r0SiO = 0;
	
	def sio2DensityToCellSize(this, rho):
		rho0 = 0.60224;
		rho = rho*rho0;
		m = 16*this.mO + 8*this.mSi;
		return (m/rho)**(1/3);
	#end
	
	def h2oDensityToCellSize(this, rho):
		rho0 = 0.60224;
		rho = rho*rho0;
		m = 2*this.mH + this.mO;
		return (m/rho)**(1/3);
	#end
#end