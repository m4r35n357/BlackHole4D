/**
 * 
 */
package uk.me.doitto;

import static uk.me.doitto.Integrator.STORMER_VERLET_10;
import static uk.me.doitto.Integrator.STORMER_VERLET_2;
import static uk.me.doitto.Integrator.STORMER_VERLET_4;
import static uk.me.doitto.Integrator.STORMER_VERLET_6;
import static uk.me.doitto.Integrator.STORMER_VERLET_8;

/**
 * @author ian
 * Particle trajectories in the Kerr spacetime and Boyer-Lindquist coordinates
 */
public final class KerrMotion {
	
	private static final double TWOPI = 2.0 * Math.PI;
	
	private final double M, a, horizon, mu2, E, L, CC, time, step, a2, lmae2; // constants for this spacetime
	
	private double r2, ra2, sth, cth, sth2, cth2, sth3, cth3, csth, sigma, delta, R, P1, P2, THETA, TH;  // intermediate variables
	
	private double tau, t, r, theta, phi, rDot, thDot, x, y, z; // coordinates etc.
	
	private Integrator symplectic;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double mass, double spin, double m, double E, double L, double C, double r, double th, double ph, double T, double ts, int order) {
		M = mass;
		a = spin;
		mu2 = m;
		this.E = E;
		this.L = L;
		CC = C;
		this.r = r;
		theta = th;
		phi = ph;
		time = T;
		step = ts;
		switch (order) {
			case 2: symplectic = STORMER_VERLET_2; break;
			case 4: symplectic = STORMER_VERLET_4; break;
			case 6: symplectic = STORMER_VERLET_6; break;
			case 8: symplectic = STORMER_VERLET_8; break;
			case 10: symplectic = STORMER_VERLET_10; break;
		}
		a2 = a * a;
		horizon = M * (1.0 + Math.sqrt(1.0 - a2));
		lmae2 = (L - a * E) * (L - a * E);
	}

	private void updateIntermediates () {
		r2 = r * r;
		ra2 = r2 + a2;
		double limit = 0.0;
		double tmp = Math.sin(theta);
		sth = Math.abs(tmp) > limit ? tmp : limit * Math.signum(tmp);
		cth = Math.cos(theta);
		sth2 = sth * sth;
		assert sth2 > 0.0 : "ZERO DIVISOR: sin(theta), theta = " + theta;
		sth3 = sth2 * sth;
		cth2 = cth * cth;
		cth3 = cth2 * cth;
		csth = cth * sth;
		sigma = r2 + a2 * cth2;
		assert sigma > 0.0 : "ZERO DIVISOR: sigma, r = " + r + ", theta = " + theta;
		delta = ra2 - 2.0 * M * r;
		assert delta > 0.0 : "ZERO DIVISOR: delta, r = " + r;
		P1 = ra2 * E - a * L;  // MTW eq.33.33b
		P2 = mu2 * r2 + lmae2 + CC;
		R = P1 * P1 - delta * P2;
		assert R > 0.0 : "R potential = " + R;
		TH = a2 * (mu2 - E * E) + L * L /sth2;
		THETA = CC - cth2 * TH;
		assert THETA > 0.0 : "THETA potential = " + THETA;
	}
	
	private double uT () {  // MTW eq.33.32d
		return ra2 * P1 / delta - a * (a * E * sth2 - L);
	}
	
	private double uPh () {  // MTW eq.33.32c
		return a * P1 / delta - (a * E - L / sth2);
	}
	
	private double hR () {
		return 10.0 * Math.log10(Math.abs(rDot * rDot - R) / 2.0);
	}
	
	private double hTh () {
		return 10.0 * Math.log10(Math.abs(thDot * thDot - THETA) / 2.0);
	}
	
	private double h () {
		return 10.0 * Math.log10(Math.abs(rDot * rDot - R) / 2.0 + Math.abs(thDot * thDot - THETA) / 2.0);
	}
	
	void updateQ (double c) {
		double cStep = c * step;
		t -= cStep * uT();
		r += cStep * rDot;
		theta += cStep * thDot;
		phi += cStep * uPh();
		updateIntermediates();
	}
	
	void updateP (double c) {
		double cStep = c * step;
		rDot += cStep * (2.0 * r * E * P1 - P2 * (r - M) - mu2 * r * delta);
		thDot += cStep * (csth * TH + L * L * cth3 / sth3);
	}
	
	public void simulate () {
		updateIntermediates();
		rDot = (R >= 0.0) ? Math.sqrt(R) : -Math.sqrt(-R);  // MTW eq.33.32b and 33.33c
		thDot = (THETA >= 0.0) ? Math.sqrt(THETA) : -Math.sqrt(-THETA);  // MTW eq.33.32a and 33.33a
		symplectic.init();
		do {
			double ra = Math.sqrt(ra2);
			System.out.printf("{\"mino\":%.9e, \"tau\":%.9e, \"H\":%.1f, \"HR\":%.1f, \"HTh\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"uT\":%.9e, \"uR\":%.9e, \"uTh\":%.9e, \"uPh\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", tau, sigma * tau, h(), hR(), hTh(), - t, r, theta, phi, uT(), rDot, thDot, uPh(), ra * sth * Math.cos(phi), ra * sth * Math.sin(phi), r * cth);
			tau += step;
			symplectic.solve(this);
		} while (r > horizon && tau <= time);  // outside horizon and in proper time range
	}
}
