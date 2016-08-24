class Parameters{
	public:
		void parameters();
		double twoBodyForce(Particle *p1, Particle *p2);
		double twoBodyOxygenForce(Particle *p1, Particle *p2);
		double threeBodyForce(Particle *p1, Particle *p2, Particle *p3);
	private:
		void Tabulate2BPotential(int type1, int type2, double *f, double *V);
		double rc, rc2;
		double types, typeSi, typeOs, typeH, typeOh;
		double L0, E0, m0, v0, t0, T0, q0;
		double *Z, *m, *Rb, *Db;
		double **H, **D, **W, **eta, **r1si, **r4si, **r0;
		double ***B, ***c0, ***xi;
};