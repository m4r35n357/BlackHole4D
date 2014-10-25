/*
 * Copyright Ian Smith 2014
 */
package uk.me.doitto;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealVector;

/**
 * @author ian
 *
 */
public class InitialConditions {
	
	private enum Trajectory {
		PARTICLE, LIGHT;
	}
	
	private enum Spin {
		RETROGRADE, ZERO, PROGRADE;
	}
	
	private final double M = 1.0, mu, a, r0, r1, th0, factorL, time = 20.0, step = 0.001, tolerance = 1.0e-6;
	
	private double E = 1.0, L = 2.0, Q = 0.0;
	
	private int order;
	
	public InitialConditions(Trajectory t, double rMin, double rMax, double thetaMin, Spin a, double factorL, Integrator i) {
		this.mu = (t == Trajectory.PARTICLE) ? 1.0 : 0.0;
		boolean singular = abs(rMax - rMin) > 2.0 * tolerance;
		this.r0 = singular ? rMin: rMin - tolerance;
		this.r1 = singular ? rMax: rMax + tolerance;
		this.th0 = thetaMin > 0.01 ? thetaMin:  0.01;
		switch (a) {
			case RETROGRADE: this.a = -1.0; break;
			case ZERO: this.a = 0.0; break;
			case PROGRADE: this.a = 1.0; break;
		default: this.a = 0.0; break;
		}
		this.factorL = factorL;
		switch (i) {
			case STORMER_VERLET_2: this.order = 2; break;
			case STORMER_VERLET_4: this.order = 4; break;
			case STORMER_VERLET_6: this.order = 6; break;
			case STORMER_VERLET_8: this.order = 8; break;
			case STORMER_VERLET_10: this.order = 10; break;
			default: this.order = 2; break;
		}
	}

	private double rDot (double r) {
		return ((r * r + a * a) * E - a * L) * ((r * r + a * a) * E - a * L) - (r * r - 2.0 * M * r + a * a) * (mu * mu * r * r + (L - a * E) * (L - a * E) + Q);
	}
	
	private double thDot (double theta) {
		return Q - cos(theta) * cos(theta) * (a * a * (mu * mu - E * E) + L * L / (sin(theta) * sin(theta)));
	}
	
	private RealVector qDot () {
		return new ArrayRealVector(new double[] { rDot(r0), rDot(r1), thDot(th0) }, false);
	}
	
	private Array2DRowRealMatrix jacobian () {
		double p0 = E * (r0 * r0 + a * a) - a * L;
		double p1 = E * (r1 * r1 + a * a) - a * L;
		double delta0 = r0 * r0 - 2.0 * M * r0 + a * a;
		double delta1 = r1 * r1 - 2.0 * M * r1 + a * a;
		double l_ae = (L - a * E);
		return new Array2DRowRealMatrix(new double[][] {
			{ 2.0 * (r0 * r0 + a * a) * p0 + 2.0 * a * l_ae * delta0, - 2.0 * a * p0 - 2.0 * l_ae * delta0, - delta0 },
			{ 2.0 * (r1 * r1 + a * a) * p1 + 2.0 * a * l_ae * delta1, - 2.0 * a * p1 - 2.0 * l_ae * delta1, - delta1 }, 
			{ 2.0 * cos(th0) * cos(th0) * a * a * E, - 2.0 * cos(th0) * cos(th0) * L / (sin(th0) * sin(th0)), 1.0 } }, false);
	}
	
	private void solve () {
		while (qDot().dotProduct(qDot()) > 1.0e-18) {
			RealVector correction = new LUDecomposition(jacobian()).getSolver().solve(qDot());
			E -= correction.getEntry(0);
			L -= correction.getEntry(1);
			Q -= correction.getEntry(2);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		InitialConditions ic = new InitialConditions(Trajectory.PARTICLE, 12.0, 12.0, 0.0, Spin.PROGRADE, 1.0, Integrator.STORMER_VERLET_10);
		ic.solve();
		new KerrMotion(ic.M, ic.a, ic.mu, ic.E, ic.L * ic.factorL, ic.Q, ic.r1, ic.th0, ic.time, ic.step, ic.order).simulate();
		System.out.println("");
		System.out.println("{ \"M\" : " + ic.M + ",");
		System.out.println("  \"a\" : " + ic.a + ",");
		System.out.println("  \"mu\" : " + ic.mu + ",");
		System.out.println("  \"E\" : " + ic.E + ",");
		System.out.println("  \"Lz\" : " + ic.L * ic.factorL + ",");
		System.out.println("  \"C\" : " + ic.Q + ",");
		System.out.println("  \"r\" : " + ic.r1 + ",");
		System.out.println("  \"theta\" : " + ic.th0 + ",");
		System.out.println("  \"time\" : " + ic.time + ",");
		System.out.println("  \"step\" : " + ic.step + ",");
		System.out.println("  \"integratorOrder\" : " + ic.order);
		System.out.println("}");
	}
}
