/**
 * 
 */
package uk.me.doitto;

import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.log10;
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
	
	private final double M, a, horizon, mu2, E, L, CC, duration, ts, a2, lmae2; // constants for this spacetime
	
	private double r2, ra2, sth, cth, sth2, cth2, delta, R, P1, P2, THETA, TH;  // intermediate variables
	
	private double mino, t, r, th, ph, rDot, thDot, eCum, e, eR, eTh; // coordinates etc.
	
	private Integrator symplectic;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double mass, double spin, double m, double E, double L, double C, double r, double th, double T, double ts, int order) {
		M = mass;
		a = spin;
		mu2 = m;
		this.E = E;
		this.L = L;
		CC = C;
		this.r = r;
		this.th = th;
		duration = T;
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
		sth = sin(th);
		cth = cos(th);
		sth2 = sth * sth;
		cth2 = cth * cth;
		delta = ra2 - 2.0 * M * r;
		P1 = ra2 * E - a * L;  // MTW eq.33.33b
		P2 = mu2 * r2 + lmae2 + CC;
		R = P1 * P1 - delta * P2;
		TH = a2 * (mu2 - E * E) + L * L /sth2;
		THETA = CC - cth2 * TH;
	}
	
	private void errors () {
		double e_r = abs(rDot * rDot - R) / 2.0;
		double e_th = abs(thDot * thDot - THETA) / 2.0;
		eR = 10.0 * log10(e_r + 1.0e-18);
		eTh = 10.0 * log10(e_th + 1.0e-18);
		e =  10.0 * log10(e_r + e_th + 1.0e-18);
		eCum += e_r + e_th;
	}
	
	private void update_T_Phi () {
		t -= ts * (ra2 * P1 / delta - a * (a * E * sth2 - L));  // MTW eq.33.32d
		ph += ts * (a * P1 / delta - a * E + L / sth2);  // MTW eq.33.32c
	}
	
	void updateQ (double c) {
		r += c * ts * rDot;
		th += c * ts * thDot;
		updateIntermediates();
	}
	
	void updateP (double c) {
		rDot += c * ts * (2.0 * r * E * P1 - P2 * (r - M) - mu2 * r * delta);  // see Maxima file bh.wxm, Mino Time
		thDot += c * ts * (cth * sth * TH + L * L * cth2 * cth / (sth2 * sth));  // see Maxima file bh.wxm, Mino Time
	}
	
	public double simulate () {
		updateIntermediates();
		rDot = -sqrt(R >= 0.0 ? R: 0.0);  // MTW eq.33.32b and 33.33c
		thDot = -sqrt(THETA >= 0.0 ? THETA: 0.0);  // MTW eq.33.32a and 33.33a
		symplectic.init();
		do {
			errors();
			double ra = sqrt(ra2);
			System.out.printf("{\"mino\":%.9e, \"tau\":%.9e, \"E\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"R\":%.9e, \"THETA\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n",
					mino, (r2 + a2 * cth2) * mino, e, eR, eTh, - t, r, th, ph, R, THETA, ra * sth * cos(ph), ra * sth * sin(ph), r * cth);
			mino += ts;
			update_T_Phi();
			symplectic.solve(this);
		} while (r > horizon && mino <= duration);  // outside horizon and in proper time range
		return eCum;
	}
}
