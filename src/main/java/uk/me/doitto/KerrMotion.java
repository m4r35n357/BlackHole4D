/**
 * 
 */
package uk.me.doitto;

import static uk.me.doitto.Integrator.STORMER_VERLET_10;
import static uk.me.doitto.Integrator.STORMER_VERLET_2;
import static uk.me.doitto.Integrator.STORMER_VERLET_4;
import static uk.me.doitto.Integrator.STORMER_VERLET_6;
import static uk.me.doitto.Integrator.STORMER_VERLET_8;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.json.simple.JSONObject;
import org.json.simple.JSONValue;

/**
 * @author ian
 * Geodesics in the Kerr spacetime and Boyer-Lindquist coordinates
 */
public class KerrMotion {
	
	static final double TWOPI = 2.0 * Math.PI;
	
	static final double mu = 1.0;
	
	static final double mu2 = mu * mu;
	
	final double M, a, E, Lz, C, step, a2, f1, f12;
	
	private double r2, ra2, ra, sth, cth, sth2, cth2, sth3, cth3, csth, sigma, sigma2, sigma3, delta,  P, R, THETA, f2;
	
	double tau, t, r, theta, phi, rDot, thetaDot, x, y, z;
	
	private final Integrator integrator;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double mass, double spin, double E, double L, double C, double t, double r, double theta, double phi, double step, int order) {
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
		this.f12 = this.f1 * this.f1;
		updateIntermediates(r, theta);
		this.rDot = uR();
		this.thetaDot = uTh();
		switch (order) {
		case 2:
			integrator = STORMER_VERLET_2;
			break;
		case 4:
			integrator = STORMER_VERLET_4;
			break;
		case 6:
			integrator = STORMER_VERLET_6;
			break;
		case 8:
			integrator = STORMER_VERLET_8;
			break;
		case 10:
			integrator = STORMER_VERLET_10;
			break;
		default:
			integrator = STORMER_VERLET_4;
			break;
		}
		integrator.init();
	}

	private void updateIntermediates (double r, double theta) {
		r2 = r * r;
		ra2 = r2 + a2;
		ra = Math.sqrt(ra2);
		sth = Math.sin(theta);
		cth = Math.cos(theta);
		sth2 = sth * sth;
		sth3 = sth2 * sth;
		cth2 = cth * cth;
		cth3 = cth2 * cth;
		csth = cth * sth;
		sigma = r2 + a2 * cth2;
		assert sigma > 0.0 : "delta = " + sigma;
		sigma2 = sigma * sigma;
		sigma3 = sigma2 * sigma;
		delta = ra2 - 2.0 * M * r;
		assert delta > 0.0 : "delta = " + delta;
		P = ra2 * E - a * Lz;  // MTW eq.33.33b
		f2 = a2 * (mu2 - E * E) + Lz * Lz /sth2;
		R = P * P - delta * (mu2 * r2 + f1 * f1 + C);
		THETA = C - cth2 * f2;
	}
	
	private double correctTheta (double theta) {
		double newTheta = theta;
		if (theta < 0.0) {
			newTheta = - theta;
		} else if (theta > Math.PI) {
			newTheta = Math.PI - (theta % TWOPI);
		}
		return newTheta;
	}
	
	public boolean outsideHorizon () {
		return r > M * (1.0 + Math.sqrt(1.0 - a * a));
	}
	
	private double uT () {  // MTW eq.33.32d
		return (ra2 * P / delta - a * (a * E * sth2 - Lz)) / sigma;
	}
	
	private double uR () {  // MTW eq.33.32b and 33.33c
		double R = Math.abs(this.R);
		assert R >= 0.0 : "NEGATIVE SQUARE ROOT: R = " + R;
		return Math.sqrt(R) / sigma;
	}
	
	private double uTh () {  // MTW eq.33.32a and 33.33a
		double THETA = Math.abs(this.THETA);
		assert THETA >= 0.0 : "NEGATIVE SQUARE ROOT: THETA = " + THETA;
		return Math.sqrt(THETA) / sigma;
	}
	
	private double uPh () {  // MTW eq.33.32c
		return (a * P / delta - (a * E - Lz / sth2)) / sigma;
	}
	
	void updateQ (double c) {
		double tmp = c * step / mu;
		r += rDot * tmp;
		theta = correctTheta(theta + thetaDot * tmp);
		updateIntermediates(r, theta);
	}
	
	void updateP (double c) {
		double tmp = c * step;
		rDot += tmp * ((2.0 * r * E * P - mu2 * r * delta - (f12 + C + mu2 * r2) * (r - M)) / sigma2 - 2.0 * r * (R + THETA) / sigma3);
		thetaDot += tmp * ((csth * f2 + Lz * Lz * cth3 / sth3) / sigma2 + 2.0 * csth * a2 * (R + THETA) / sigma3);
	}
	
	public double v4n () {
		double h1 = (uT() - a * sth2 * uPh());
		double h2 = (ra2 * uPh() - a * uT());
		return - delta / sigma * h1 * h1 + sth2 / sigma * h2 * h2 + sigma / delta * uR() * uR() + sigma * uTh() * uTh();  // based on MTW eq. 33.2
	}
	
	public double hamiltonian () {
		return 0.5 * (rDot * rDot + thetaDot * thetaDot - (R + THETA) / sigma2);
	}
	
	public double iterateSymplectic () {
		integrator.solve(this);
		tau += step;
		t += uT() * step;
		phi = (phi + uPh() * step) % TWOPI;
		x = ra * sth * Math.cos(phi);
		y = ra * sth * Math.sin(phi);
		z = r * cth;
		return v4n();
	}
	
	/**
	 * Read initial conditions from a JSON-formatted file
	 * @param fileName the path to the file
	 * @return a Symplectic instance
	 */
	public static KerrMotion icJson (String fileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(fileName)));
		String data = "";
		String line = bufferedReader.readLine();
		while (line != null) {
			data += line;
			line = bufferedReader.readLine();
		}
		bufferedReader.close();
		JSONObject ic = (JSONObject)JSONValue.parse(data);
		return new KerrMotion ((Double)ic.get("M"), (Double)ic.get("a"), (Double)ic.get("E"), (Double)ic.get("Lz"), (Double)ic.get("C"), (Double)ic.get("t"), (Double)ic.get("r"), (Double)ic.get("theta"), (Double)ic.get("phi"), (Double)ic.get("step"), ((Long)ic.get("integratorOrder")).intValue());
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main (String[] args) throws IOException {
		KerrMotion st;
		if (args.length == 1) {
			st = KerrMotion.icJson(args[0]);
		} else {
			System.err.println("Missing file name, giving up!");
			return;
		}
//		KerrMotion st = new KerrMotion(1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 1.0 / 16.0);
//		KerrMotion st = new KerrMotion(1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 1.0 / 16.0);  // GOOD, don't touch!
//		KerrMotion st = new KerrMotion(1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 3.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 1.0 / 16.0);
//		KerrMotion st = new KerrMotion(1.0, 1.0, 0.989352727272727, -4.683, 0.0, 0.0, 12.201, Math.PI / 2.0, 0.0, 1.0 / 32.0);
//		KerrMotion st = new KerrMotion(1.0, 0.0, 1.0, 4.0, 0.0, 0.0, 4.0, Math.PI / 2.0, 0.0, 1.0 / 4.0);
//		KerrMotion st = new KerrMotion(1.0, 1.0, 0.966, 4.066, 2.0, 0.0, 17.488, Math.PI / 2.0, 0.0, 1.0 / 16.0);
		double v4Norm;
		double h;
		while (st.outsideHorizon()) {
			v4Norm = st.iterateSymplectic();
			h = st.hamiltonian();
			System.out.printf("{\"V2\":%.9e, \"H\":%.1f, \"tau\":%.9e, \"t\":%.9e, \"r\":%.9e, \"theta\":%.9e, \"phi\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", -v4Norm, 10.0 * Math.log10(Math.abs(h)), st.tau, st.t, st.r, st.theta, st.phi, st.x, st.y, st.z);
//			System.out.printf("{\"V2\":%.9e, \"H\":%.9e, \"tau\":%.9e, \"t\":%.9e, \"r\":%.9e, \"theta\":%.9e, \"phi\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", -v4Norm, h, st.tau, st.t, st.r, st.theta, st.phi, st.x, st.y, st.z);
		}
	}
}
