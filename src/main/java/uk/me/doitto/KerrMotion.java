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
 * Particle trajectories in the Kerr spacetime and Boyer-Lindquist coordinates
 */
public class KerrMotion {
	
	private static final double TWOPI = 2.0 * Math.PI;
	
	private final double M, a, horizon, mu, mu2, E, Lz, CC, time, step, a2, f1, f12; // constants for this spacetime
	
	private double r2, ra2, ra, sth, cth, sth2, cth2, sth3, cth3, csth, sigma, sigma2, sigma3, delta, P, R, THETA, R_THETA, f2;  // intermediate variables
	
	private double tau, t, r, theta, phi, rDot, thetaDot, x, y, z; // coordinates etc.
	
	private final Integrator symplectic;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double mass, double spin, double m, double E, double L, double C, double t, double r, double th, double ph, double T, double ts, int order) {
		M = mass;
		a = spin;
		horizon = M * (1.0 + Math.sqrt(1.0 - a * a));
		mu = m;
		mu2 = mu * mu;
		this.E = E;
		Lz = L;
		CC = C;
		this.t = t;
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
		R = P * P - delta * (mu2 * r2 + f1 * f1 + CC);
		f2 = a2 * (mu2 - E * E) + Lz * Lz /sth2;
		THETA = CC - cth2 * f2;
		R_THETA = R + THETA;
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
	
	private double correctTheta (double theta) {
		double newTheta = theta;
		if (theta < 0.0) {
			newTheta = - theta;
		} else if (theta > Math.PI) {
			newTheta = Math.PI - (theta % TWOPI);
		}
		return newTheta;
	}
	
	private double v2 () {
		double h1 = (uT() - a * sth2 * uPh());
		double h2 = (ra2 * uPh() - a * uT());
		return - delta / sigma * h1 * h1 + sth2 / sigma * h2 * h2 + sigma / delta * uR2() + sigma * uTh2();  // based on MTW eq. 33.2
	}
	
	private double pH () {
		return 10.0 * Math.log10(Math.abs(0.5 * (rDot * rDot + thetaDot * thetaDot - R_THETA / sigma2)));
	}
	
	void updateQ (double c) {
		double tmp = c * step / mu;
		t += uT() * tmp;
		r += rDot * tmp;
		theta = correctTheta(theta + thetaDot * tmp);
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
		rDot = Math.sqrt(uR2()) >= 0.0 ? Math.sqrt(uR2()) : 0.0;
		thetaDot = Math.sqrt(uTh2()) >= 0.0 ? Math.sqrt(uTh2()) : 0.0;
		symplectic.init();
		do {
			x = ra * sth * Math.cos(phi);
			y = ra * sth * Math.sin(phi);
			z = r * cth;
			System.out.printf("{\"v2\":%.3f, \"H\":%.1f, \"tau\":%.9e, \"t\":%.9e, \"r\":%.9e, \"theta\":%.9e, \"phi\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n", - v2(), pH(), - tau, - t, r, theta, phi, x, y, z);
			tau += step;
			symplectic.solve(this);
		} while (r > horizon && - tau <= time);  // outside horizon and in proper time range
	}
	
	/**
	 * Read initial conditions from a JSON-formatted parameter file using Google's SimpleJSON library
	 * @param fileName the path to the file
	 * @return a KerrMotion instance
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
