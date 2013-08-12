/**
 * 
 */
package uk.me.doitto;
/**
 * @author ian
 *
 */
public class KerrMotion {
	
	static final double TWOPI = 2.0 * Math.PI;
	
	static final double c = 299792458.0;
	
	static final double G = 6.67398e-11;
	
	static final double mu2 = 1.0;
	
	final double M, a, E, Lz, C, step, a2, f1;
	
	private double r2, ra2, ra, sth, cth, sth2, cth2, sigma, delta, P, sgn_r, sgn_r_set_value, sgn_theta, sgn_theta_set_value;
	
	private boolean sgn_r_set, sgn_theta_set;

	double tau, t, r, theta, phi, tDot, rDot, thetaDot, phiDot, x, y, z;
	
	/**
	 * 
	 */
	public KerrMotion(double mass, double spin, double E, double L, double C, double t, double r, double theta, double phi, double step) {
		this.M = mass;
		this.a = spin;
		this.E = E;
		this.Lz = L;
		this.C = C;
		this.tau = 0.0;
		this.t = t;
		this.r = r;
		this.theta = theta;
		this.phi = phi;
		this.step = - step;
		this.a2 = spin * spin;
		this.f1 = L - spin * E;
		updateIntermediates(r, theta);
	}

	static KerrMotion constantR (double mass, double spin, double E, double t, double r, double theta, double phi, double step) {
		double B = 1 / (spin * (mass - r));  // TODO this is a BIG problem for a = 0.0 !!!
		double a2 = spin * spin;
		double r2 = r * r;
		double r4 = r2 * r2;
		double M2r = 2 * mass * r;
		double E2 = E * E;
		double M2 = mass * mass;
		/* Calculating Phi and Q*/
		double Lz = B * (mass * E * (a2 - r2) + (a2 + r2 - M2r) * Math.sqrt(r2 * E2 + r * (mass- r)));
		double Q = r2 * B * B * ((1 - E2) * r4 + (4 * mass * E2 - 5 * mass) * r2 * r + (8 - 5 * E2) * M2 * r2 +
			(2 * a2 * E2 - a2 - 4 * M2) * mass * r + a2 * M2 + 2 * mass * E * (a2 + r2 - M2r) * Math.sqrt(r * (mass - r + r * E2)));
		return new KerrMotion(mass, spin, E, Lz, Q, t, r, theta, phi, step);
	}
	
	private void updateIntermediates (double r, double theta) {
		r2 = r * r;
		ra2 = r2 + a2;
		ra = Math.sqrt(ra2);
		sth = Math.sin(theta);
		cth = Math.cos(theta);
		sth2 = sth * sth;
		assert sth2 > 0.0 : "Sin^2(theta) = " + sth2;
		cth2 = cth * cth;
		sigma = r2 + a2 * cth2;
		assert sigma > 0.0 : "delta = " + sigma;
		delta = ra2 - 2.0 * M * r;
		assert delta > 0.0 : "delta = " + delta;
		P = E * ra2 - Lz * a;
	}
	
	public boolean outsideHorizon () {
		return r > M * (1.0 + Math.sqrt(1.0 - a * a));
	}
	
	private double uT () {
		return (ra2 * P / delta - a * (a * E * sth2 - Lz)) / sigma;
	}
	
	private double uR () {
		double R = Math.abs(P * P - delta * (mu2 * r2 + f1 * f1 + C));
		assert R >= 0.0 : "NEGATIVE SQUARE ROOT: R = " + R;
		return Math.sqrt(R) / sigma;
	}
	
	private double uTh () {
		double THETA = Math.abs(C - cth2 * (a2 * (mu2 - E * E) + Lz * Lz / sth2));
		assert THETA >= 0.0 : "NEGATIVE SQUARE ROOT: THETA = " + THETA;
		return Math.sqrt(THETA) / sigma;
	}
	
	private double uPh () {
		return (a * P / delta - (a * E - Lz / sth2)) / sigma;
	}
	
	private double updateCoordinates (double deltaT, double deltaR, double deltaTh, double deltaPh) {
		tau += step;
		t += deltaT;		
		r += deltaR;
		theta += deltaTh;
		phi += deltaPh;
//		if (Math.abs(deltaR) < 10 * step * step && !sgn_r_set) {
//			sgn_r = -sgn_r;
//			sgn_r_set = true;			/* lock to prevent further sign changes */
//			sgn_r_set_value = r;
//		}
//		if (sgn_r_set && Math.abs(sgn_r_set_value - r) > step) {
//			sgn_r_set = false;			/* unlock */
//		}
//		if (Math.abs(deltaTh) < 10 * step * step && !sgn_theta_set) {
//			sgn_theta = -sgn_theta;
//			sgn_theta_set = true;			/* lock to prevent further sgn changes */
//			sgn_theta_set_value = theta;
//		}
//		if (sgn_theta_set && Math.abs(sgn_theta_set_value - theta) > step) {
//			sgn_theta_set = false;			/* unlock */
//		}
		updateIntermediates(r, theta);
		x = ra * sth * Math.cos(phi);
		y = ra * sth * Math.sin(phi);
		z = r * cth;
		double h1 = (uT() - a * sth2 * uPh());
		double h2 = (ra2 * uPh() - a * uT());
		return - delta / sigma * h1 * h1 + sth2 / sigma * h2 * h2 + sigma / delta * uR() * uR() + sigma * uTh() * uTh();
	}
	
	public double iterateEuler () {
		return updateCoordinates(uT() * step, uR() * step, uTh() * step, uPh() * step);
	}
	
	private double rk4Up (double k1, double k2, double k3, double k4) {
		return (k1 + 2.0 * (k2 + k3) + k4) * step / 6.0;
	}
	
	public double iterateRk4 () {
		double tK1 = uT();
		double rK1 = uR();
		double thK1 = uTh();
		double phK1 = uPh();
		updateIntermediates(r + rK1 * step / 2.0, theta + thK1 * step / 2.0);  // RK4 stage 2
		double tK2 = uT();
		double rK2 = uR();
		double thK2 = uTh();
		double phK2 = uPh();
		updateIntermediates(r + rK2 * step / 2.0, theta + thK2 * step / 2.0);// RK4 stage 3
		double tK3 = uT();
		double rK3 = uR();
		double thK3 = uTh();
		double phK3 = uPh();
		updateIntermediates(r + rK3 * step, theta + thK3 * step);// RK4 stage 4
		double tK4 = uT();
		double rK4 = uR();
		double thK4 = uTh();
		double phK4 = uPh();
		return updateCoordinates(rk4Up(tK1, tK2, tK3, tK4), rk4Up(rK1, rK2, rK3, rK4), rk4Up(thK1, thK2, thK3, thK4), rk4Up(phK1, phK2, phK3, phK4));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		System.out.println(args.length);
		KerrMotion st;
//		if (args.length == 12) {
//			st = new KerrMotion(Double.parseDouble(args[0]), Double.parseDouble(args[1]), Double.parseDouble(args[2]), Double.parseDouble(args[3]), Double.parseDouble(args[4]) ,Double.parseDouble(args[5]) ,Double.parseDouble(args[6]) ,Double.parseDouble(args[7]) ,Double.parseDouble(args[8]) ,Double.parseDouble(args[9]));
//		} else {
			st = new KerrMotion(1.0, 0.0, 0.962250448649377, 4.0, 0.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 0.03125);
//		}
		double hamiltonian;
//		KerrMotion st = KerrMotion.constantR(1.0, 1.0, 0.960, 0.0, 10.0, Math.PI / 2.0, 0.0, 0.125);
		while (st.outsideHorizon()) {
			hamiltonian = st.iterateEuler();
			System.out.printf("{\"H\":%.9e, \"tau\":%.9e, \"t\":%.9e, \"r\":%.9e, \"theta\":%.9e, \"phi\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", hamiltonian, st.tau, st.t, st.r, st.theta % TWOPI, st.phi % TWOPI, st.x, st.y, st.z);
		}
	}
}
