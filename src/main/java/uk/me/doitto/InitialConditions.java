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
	
	private double M = 1.0, mu = 1.0, E = 1.0, L = 2.0, Q = 0.0, a, r0, r1, theta0, factorL = 1.0, tolerance = 1.0e-6;
	
	public InitialConditions(double rMin, double rMax, double thetaMin, double a, double factorL) {
		boolean singular = abs(rMax - rMin) > 2.0 * tolerance;
		this.r0 = singular ? rMin: rMin - tolerance;
		this.r1 = singular ? rMax: rMax + tolerance;
		this.theta0 = thetaMin > 0.01 ? thetaMin:  0.01;
		this.a = a;
		this.factorL = factorL;
	}

	private double rDot (double r) {
		return ((r * r + a * a) * E - a * L) * ((r * r + a * a) * E - a * L) - (r * r - 2.0 * M * r + a * a) * (mu * mu * r * r + (L - a * E) * (L - a * E) + Q);
	}
	
	private double thetaDot (double theta) {
		return Q - cos(theta) * cos(theta) * (a * a * (mu * mu - E * E) + L * L / (sin(theta) * sin(theta)));
	}
	
	private RealVector qDot () {
		return new ArrayRealVector(new double[] { rDot(r0), rDot(r1), thetaDot(theta0) }, false);
	}
	
	private Array2DRowRealMatrix getJacobian () {
		double p0 = E * (r0 * r0 + a * a) - a * L;
		double p1 = E * (r1 * r1 + a * a) - a * L;
		double delta0 = r0 * r0 - 2.0 * M * r0 + a * a;
		double delta1 = r1 * r1 - 2.0 * M * r1 + a * a;
		double l_ae = (L - a * E);
		return new Array2DRowRealMatrix(new double[][] {
			{ 2.0 * (r0 * r0 + a * a) * p0 + 2.0 * a * l_ae * delta0, -2.0 * a * p0 - 2.0 * l_ae * delta0, - delta0 },
			{ 2.0 * (r1 * r1 + a * a) * p1 + 2.0 * a * l_ae * delta1, -2.0 * a * p1 - 2.0 * l_ae * delta1, - delta1 }, 
			{ 2.0 * cos(theta0) * cos(theta0) * a * a * E, - 2.0 * cos(theta0) * cos(theta0) * L / (sin(theta0) * sin(theta0)), 1.0 } }, false);
	}
	
	private void solve () {
		while (qDot().dotProduct(qDot()) > 1.0e-18) {
			RealVector correction = new LUDecomposition(getJacobian()).getSolver().solve(qDot());
			E -= correction.getEntry(0);
			L -= correction.getEntry(1);
			Q -= correction.getEntry(2);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		InitialConditions ic = new InitialConditions(12.0, 12.0, 0.0, 1.0, 1.0);
		ic.solve();
		new KerrMotion(1.0, ic.a, 1.0, ic.E, ic.factorL * ic.L, ic.Q, sqrt(ic.r0 * ic.r1), PI / 2.0, 50.0, 0.001, 8).simulate();
		System.out.println("");
		System.out.println("{ \"M\" : 1.0,");
		System.out.println("  \"a\" : " + ic.a + ",");
		System.out.println("  \"mu\" : 1.0,");
		System.out.println("  \"E\" : " + ic.E + ",");
		System.out.println("  \"Lz\" : " + ic.factorL * ic.L + ",");
		System.out.println("  \"C\" : " + ic.Q + ",");
		System.out.println("  \"r\" : " + sqrt(ic.r0 * ic.r1) + ",");
		System.out.println("  \"theta\" : " + PI / 2.0 + ",");
		System.out.println("  \"time\" : 20.0,");
		System.out.println("  \"step\" : 0.001,");
		System.out.println("  \"integratorOrder\" : 8");
		System.out.println("}");
	}
}
