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
	
	private final double M, a, horizon, mu, mu2, E, Lz, CC, time, step, a2, f1, f12; // constants for this spacetime
	
	private double r2, ra2, ra, sth, cth, sth2, cth2, sth3, cth3, csth, sigma, sigma2, sigma3, delta, P, R, THETA, R_THETA, f2;  // intermediate variables
	
	private double tau, t, r, theta, phi, rDot, thetaDot, x, y, z; // coordinates etc.
	
	private final Integrator symplectic;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double mass, double spin, double m, double E, double L, double C, double r, double th, double ph, double T, double ts, int order) {
		M = mass;
		a = spin;
		mu = m;
		this.E = E;
		Lz = L;
		CC = C;
		this.r = r;
		theta = th;
		phi = ph;
		time = T;
		step = - ts;
		switch (order) {
		case 2:
			symplectic = STORMER_VERLET_2;
			break;
		case 4:
			symplectic = STORMER_VERLET_4;
			break;
		case 6:
			symplectic = STORMER_VERLET_6;
			break;
		case 8:
			symplectic = STORMER_VERLET_8;
			break;
		case 10:
			symplectic = STORMER_VERLET_10;
			break;
		default:
			symplectic = STORMER_VERLET_4;
			break;
		}
		a2 = a * a;
		horizon = M * (1.0 + Math.sqrt(1.0 - a2));
		mu2 = mu * mu;
		f1 = L - a * E;
		f12 = f1 * f1;
	}

	private void updateIntermediates (double r, double theta) {
		r2 = r * r;
		ra2 = r2 + a2;
		ra = Math.sqrt(ra2);
		sth = Math.sin(theta);
		cth = Math.cos(theta);
		sth2 = sth * sth;
		assert sth2 > 0.0 : "ZERO DIVISOR: sin(theta), theta = " + theta;
		sth3 = sth2 * sth;
		cth2 = cth * cth;
		cth3 = cth2 * cth;
		csth = cth * sth;
		sigma = r2 + a2 * cth2;
		assert sigma > 0.0 : "ZERO DIVISOR: sigma, r = " + r + ", theta = " + theta;
		sigma2 = sigma * sigma;
		sigma3 = sigma2 * sigma;
		delta = ra2 - 2.0 * M * r;
		assert delta > 0.0 : "ZERO DIVISOR: delta, r = " + r;
		P = ra2 * E - a * Lz;  // MTW eq.33.33b
		R = P * P - delta * (mu2 * r2 + f12 + CC);
		f2 = a2 * (mu2 - E * E) + Lz * Lz /sth2;
		THETA = CC - cth2 * f2;
		R_THETA = R + THETA;
	}
	
	private double uT () {  // MTW eq.33.32d
		return (ra2 * P / delta - a * (a * E * sth2 - Lz)) / sigma;
	}
	
	private double uPh () {  // MTW eq.33.32c
		return (a * P / delta - (a * E - Lz / sth2)) / sigma;
	}
	
	private double v2 () {
		double h1 = (uT() - a * sth2 * uPh());
		double h2 = (ra2 * uPh() - a * uT());
		return - delta / sigma * h1 * h1 + sth2 / sigma * h2 * h2 + R / (delta * sigma) + THETA / sigma;  // based on MTW eq. 33.2
	}
	
	private double pH () {
		return 10.0 * Math.log10(Math.abs(0.5 * (rDot * rDot + thetaDot * thetaDot - R_THETA / sigma2)));
	}
	
	void updateQ (double c) {
		double tmp = c * step / mu;
		t += uT() * tmp;
		r += rDot * tmp;
		theta = (theta + thetaDot * tmp) % TWOPI;
		phi = (phi - uPh() * tmp) % TWOPI;
		updateIntermediates(r, theta);
	}
	
	void updateP (double c) {
		double tmp = c * step;  // NB factor of 0.5 cancelled out by a factor of 2.0 below
		rDot += tmp * ((2.0 * r * E * P - mu2 * r * delta - (f12 + CC + mu2 * r2) * (r - M)) / sigma2 - 2.0 * r * R_THETA / sigma3);
		thetaDot += tmp * ((csth * f2 + Lz * Lz * cth3 / sth3) / sigma2 + 2.0 * csth * a2 * R_THETA / sigma3);
	}
	
	public void simulate () {
		updateIntermediates(r, theta);
		rDot = (R / sigma2 >= 0.0) ? Math.sqrt(R / sigma2) : - Math.sqrt(- R / sigma2);  // MTW eq.33.32b and 33.33c
		thetaDot = (THETA / sigma2 >= 0.0) ? Math.sqrt(THETA / sigma2) : - Math.sqrt(- THETA / sigma2);  // MTW eq.33.32a and 33.33a
		symplectic.init();
		do {
			x = ra * sth * Math.cos(phi);
			y = ra * sth * Math.sin(phi);
			z = r * cth;
			System.out.printf("{\"tau\":%.9e, \"H\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"uT\":%.9e, \"uR\":%.9e, \"uTh\":%.9e, \"uPh\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", - tau, pH(), - t, r, theta, phi, uT(), rDot, thetaDot, uPh(), x, y, z);
			System.err.printf("{\"tau\":%.9e, \"v2\":%.3f, \"H\":%.1f}%n", - tau, - v2(), pH());
			tau += step;
			symplectic.solve(this);
		} while (r > horizon && - tau <= time);  // outside horizon and in proper time range
	}
}
