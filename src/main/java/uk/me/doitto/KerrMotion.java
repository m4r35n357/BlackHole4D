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
	
	final double M, a, mu, mu2, E, Lz, C, time, step, a2, f1, f12;
	
	private double r2, ra2, ra, sth, cth, sth2, cth2, sth3, cth3, csth, sigma, sigma2, sigma3, delta, P, R, THETA, R_THETA, f2;
	
	double tau, t, r, theta, phi, rDot, thetaDot, x, y, z;
	
	private final Integrator symplectic;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double mass, double spin, double mu, double E, double L, double C, double t, double r, double theta, double phi, double time, double step, int order) {
		this.M = mass;
		this.a = spin;
		this.mu = mu;
		this.mu2 = mu * mu;
		this.E = E;
		this.Lz = L;
		this.C = C;
		this.tau = 0.0;
		this.t = t;
		this.r = r;
		this.theta = theta;
		this.phi = phi;
		this.time = time;
		this.step = - step;
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
		symplectic.init();
		this.a2 = spin * spin;
		this.f1 = L - spin * E;
		this.f12 = this.f1 * this.f1;
		updateIntermediates(r, theta);
		this.rDot = Math.sqrt(uR2()) >= 0.0 ? Math.sqrt(uR2()) : 0.0;
		this.thetaDot = Math.sqrt(uTh2()) >= 0.0 ? Math.sqrt(uTh2()) : 0.0;
		x = ra * sth * Math.cos(phi);
		y = ra * sth * Math.sin(phi);
		z = r * cth;
		System.out.printf("{\"v2\":%.3f, \"H\":%.1f, \"tau\":%.9e, \"t\":%.9e, \"r\":%.9e, \"theta\":%.9e, \"phi\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", -v4n(), 10.0 * Math.log10(Math.abs(hamiltonian())), tau, t, r, theta, phi, x, y, z);
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
		R = P * P - delta * (mu2 * r2 + f1 * f1 + C);
		f2 = a2 * (mu2 - E * E) + Lz * Lz /sth2;
		THETA = C - cth2 * f2;
		R_THETA = R + THETA;
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
	
	private double uT () {  // MTW eq.33.32d
		return (ra2 * P / delta - a * (a * E * sth2 - Lz)) / sigma;
	}
	
	private double uR2 () {  // MTW eq.33.32b and 33.33c
		return R / sigma2;
	}
	
	private double uTh2 () {  // MTW eq.33.32a and 33.33a
		return THETA / sigma2;
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
		rDot += tmp * ((2.0 * r * E * P - mu2 * r * delta - (f12 + C + mu2 * r2) * (r - M)) / sigma2 - 2.0 * r * R_THETA / sigma3);
		thetaDot += tmp * ((csth * f2 + Lz * Lz * cth3 / sth3) / sigma2 + 2.0 * csth * a2 * R_THETA / sigma3);
	}
	
	private double v4n () {
		double h1 = (uT() - a * sth2 * uPh());
		double h2 = (ra2 * uPh() - a * uT());
		return - delta / sigma * h1 * h1 + sth2 / sigma * h2 * h2 + sigma / delta * uR2() + sigma * uTh2();  // based on MTW eq. 33.2
	}
	
	private double hamiltonian () {
		return 0.5 * (rDot * rDot + thetaDot * thetaDot - R_THETA / sigma2);
	}
	
	public void simulate () {
		while (r > M * (1.0 + Math.sqrt(1.0 - a * a)) && -tau < time) {  // outside horizon and in proper time range
			tau += step;
			symplectic.solve(this);  // symplectic integrator for r and theta
			t += uT() * step;  // euler for t
			phi = (phi + uPh() * step) % TWOPI;  // euler for phi
			x = ra * sth * Math.cos(phi);
			y = ra * sth * Math.sin(phi);
			z = r * cth;
			System.out.printf("{\"v2\":%.3f, \"H\":%.1f, \"tau\":%.9e, \"t\":%.9e, \"r\":%.9e, \"theta\":%.9e, \"phi\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", -v4n(), 10.0 * Math.log10(Math.abs(hamiltonian())), tau, t, r, theta, phi, x, y, z);
		}
	}
	
	/**
	 * Read initial conditions from a JSON-formatted file using Google's SimpleJSON library
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
		return new KerrMotion ((Double)ic.get("M"), (Double)ic.get("a"), (Double)ic.get("mu"), (Double)ic.get("E"), (Double)ic.get("Lz"), (Double)ic.get("C"), (Double)ic.get("t"), (Double)ic.get("r"), (Double)ic.get("theta"), (Double)ic.get("phi"), (Double)ic.get("time"), (Double)ic.get("step"), ((Long)ic.get("integratorOrder")).intValue());
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main (String[] args) throws IOException {
		if (args.length == 1) {
			KerrMotion.icJson(args[0]).simulate();
		} else {
			System.err.println("Missing file name, giving up!");
		}
	}
}
