/**
 * 
 */
package uk.me.doitto;

import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.log10;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
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
	
	private final double M, a, horizon, mu2, E, L, CC, time, ts, a2, lmae2; // constants for this spacetime
	
	private double r2, ra2, sth, cth, sth2, cth2, sth3, cth3, csth, sigma, delta, R, P1, P2, THETA, TH;  // intermediate variables
	
	private double tau, t, r, theta, phi, rDot, thDot, eCum, e, eR, eTh; // coordinates etc.
	
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
		this.ts = ts;
		switch (order) {
			case 2: symplectic = STORMER_VERLET_2; break;
			case 4: symplectic = STORMER_VERLET_4; break;
			case 6: symplectic = STORMER_VERLET_6; break;
			case 8: symplectic = STORMER_VERLET_8; break;
			case 10: symplectic = STORMER_VERLET_10; break;
		}
		a2 = a * a;
		horizon = M * (1.0 + sqrt(1.0 - a2));
		lmae2 = (L - a * E) * (L - a * E);
	}

	private void updateIntermediates () {
		r2 = r * r;
		ra2 = r2 + a2;
		double limit = 0.0;
		double tmp = sin(theta);
		sth = abs(tmp) > limit ? tmp : limit * signum(tmp);
		cth = cos(theta);
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
		TH = a2 * (mu2 - E * E) + L * L /sth2;
		THETA = CC - cth2 * TH;
	}
	
	private void errors () {
		double er = abs(rDot * rDot - R) / 2.0;
		double eth = abs(thDot * thDot - THETA) / 2.0;
		eR = 10.0 * log10(er + 1.0e-18);
		eTh = 10.0 * log10(eth + 1.0e-18);
		e =  10.0 * log10(er + eth + 1.0e-18);
		eCum += er + eth;
	}
	
	void update_T_Phi () {
		t -= ts * (ra2 * P1 / delta - a * (a * E * sth2 - L));  // MTW eq.33.32d
		phi += ts * (a * P1 / delta - a * E + L / sth2);  // MTW eq.33.32c
	}
	
	void updateQ (double c) {
		r += c * ts * rDot;
		theta += c * ts * thDot;
		updateIntermediates();
	}
	
	void updateP (double c) {
		rDot += c * ts * (2.0 * r * E * P1 - P2 * (r - M) - mu2 * r * delta);  // see Maxima file bh.wxm, Mino Time
		thDot += c * ts * (csth * TH + L * L * cth3 / sth3);  // see Maxima file bh.wxm, Mino Time
	}
	
	public double simulate () {
		updateIntermediates();
		rDot = -sqrt(R >= 0.0 ? R: 0.0);  // MTW eq.33.32b and 33.33c
		thDot = -sqrt(THETA >= 0.0 ? THETA: 0.0);  // MTW eq.33.32a and 33.33a
		symplectic.init();
		do {
			errors();
			double ra = sqrt(ra2);
			System.out.printf("{\"mino\":%.9e, \"tau\":%.9e, \"H\":%.1f, \"HR\":%.1f, \"HTh\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"R\":%.9e, \"THETA\":%.9e, \"uR\":%.9e, \"uTh\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", tau, sigma * tau, e, eR, eTh, - t, r, theta, phi, R, THETA, rDot, thDot, ra * sth * cos(phi), ra * sth * sin(phi), r * cth);
			tau += ts;
			update_T_Phi();
			symplectic.solve(this);
		} while (r > horizon && tau <= time);  // outside horizon and in proper time range
		return eCum;
	}
}
