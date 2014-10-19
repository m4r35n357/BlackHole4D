package uk.me.doitto;

/**
 * @author ian
 * <p>
 * Symplectic integrators of various orders and types
 */
public enum Integrator {
	/**
	 * Stormer-Verlet 2nd-order
	 */
	STORMER_VERLET_2 {
		@Override
		protected void init() {
			double[] coefficients = {1.0};
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 4th-order
	 */
	STORMER_VERLET_4 {
		@Override
		protected void init() {
			double CUBE_ROOT_2 = Math.pow(2.0, 1.0 / 3.0);
			double y = 1.0 / (2.0 - CUBE_ROOT_2);
			double[] coefficients = {y, - y * CUBE_ROOT_2};
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 6th-order
	 */
	STORMER_VERLET_6 {
		@Override
		protected void init() {
			double[] coefficients = {0.78451361047755726381949763,
									0.23557321335935813368479318,
									-1.17767998417887100694641568,
									1.31518632068391121888424973};
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 8th-order
	 */
	STORMER_VERLET_8 {
		@Override
		protected void init() {
			double[] coefficients = {0.74167036435061295344822780,
									-0.40910082580003159399730010,
									0.19075471029623837995387626,
									-0.57386247111608226665638773,
									0.29906418130365592384446354,
									0.33462491824529818378495798,
									0.31529309239676659663205666,
									-0.79688793935291635401978884};
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 10th-order
	 */
	STORMER_VERLET_10 {
		@Override
		protected void init() {
			double[] coefficients = {0.09040619368607278492161150,
									0.53591815953030120213784983,
									0.35123257547493978187517736,
									-0.31116802097815835426086544,
									-0.52556314194263510431065549,
									0.14447909410225247647345695,
									0.02983588609748235818064083,
									0.17786179923739805133592238,
									0.09826906939341637652532377,
									0.46179986210411860873242126,
									-0.33377845599881851314531820,
									0.07095684836524793621031152,
									0.23666960070126868771909819,
									-0.49725977950660985445028388,
									-0.30399616617237257346546356,
									0.05246957188100069574521612,
									0.44373380805019087955111365};
			this.coefficients = coefficients;
		}
	};

	protected double[] coefficients;
	
	/**
	 * Basic 2nd-order Stormer-Verlet step which is composed into higher order methods
	 * @param s the KerrMotion instance
	 * @param y composition coefficient
	 */
	protected final void sympBase (KerrMotion s, double y) {
		double halfY = 0.5 * y;
		s.updateQ(halfY);
		s.updateP(y);
		s.updateQ(halfY);
	}
	
	/**
	 * Set up the integrator composition coefficients
	 */
	protected abstract void init ();
	
	/**
	 * Perform one iteration step for the configured integrator
	 * @param s the KerrMotion object reference, for passing through to the Q & P update methods
	 */
	void solve (KerrMotion s) {
		int tmp = coefficients.length - 1;
		for (int i = 0; i < tmp; i++) {
			sympBase(s, coefficients[i]);
		}
		for (int i = tmp; i >= 0; i--) {
			sympBase(s, coefficients[i]);
		}
	}
}
