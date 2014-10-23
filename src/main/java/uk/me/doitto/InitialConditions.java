/**
 * 
 */
package uk.me.doitto;

import static java.lang.Math.PI;
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
	
	private double M = 1.0, mu = 1.0, E = 1.0, L = 2.0, Q = 0.0, a, rMin, rMax, thetaMin, tolerance = 1.0e-9;
	
	/**
	 * 
	 */
	public InitialConditions(double rMin, double rMax, double thetaMin, double a) {
		this.rMin = rMin - tolerance;
		this.rMax = rMax + tolerance;
		this.thetaMin = thetaMin;
		this.a = a;
	}

	private double rDot (double r) {
		return ((r * r + a * a) * E - a * L) * ((r * r + a * a) * E - a * L) - (r * r - 2.0 * M * r + a * a) * (mu * mu * r * r + (L - a * E) * (L - a * E) + Q);
	}
	
	private double thetaDot (double theta) {
		return Q - cos(theta) * cos(theta) * (a * a * (mu * mu - E * E) + L * L / (sin(theta) * sin(theta)));
	}
	
	private RealVector qDot () {
		return new ArrayRealVector(new double[] { rDot(rMin), rDot(rMax), thetaDot(thetaMin) }, false);
	}
	
	private Array2DRowRealMatrix getJacobian () {
		return new Array2DRowRealMatrix(new double[][] {
			{ 2.0 * (rMin * rMin + a * a) * (E * (rMin * rMin + a * a) - a * L) + 2.0 * a * (L - a * E) * (rMin * rMin - 2.0 * M * rMin * rMin + a * a),
				-2.0 * a * (E * (rMin * rMin + a * a) - a * L) - 2.0 * (L - a * E) * (rMin * rMin - 2.0 * M * rMin + a * a),
				- rMin * rMin + 2.0 * M * rMin - a * a },
			{ 2.0 * (rMax * rMax + a * a) * (E * (rMax * rMax + a * a) - a * L) + 2.0 * a * (L - a * E) * (rMax * rMax - 2.0 * M * rMax * rMax + a * a),
				-2.0 * a * (E * (rMax * rMax + a * a) - a * L) - 2.0 * (L - a * E) * (rMax * rMax - 2.0 * M * rMax + a * a),
				- rMax * rMax + 2.0 * M * rMax - a * a }, 
			{ 2.0 * cos(thetaMin) * cos(thetaMin) * a * a * E, - 2.0 * cos(thetaMin) * cos(thetaMin) * L / (sin(thetaMin) * sin(thetaMin)), 1.0 } }, false);
	}
	
	private void solve () {
		while (qDot().dotProduct(qDot()) > 1.0e-18) {
			RealVector correction = new LUDecomposition(getJacobian()).getSolver().solve(qDot());
			E -= correction.getEntry(0);
			L -= correction.getEntry(1);
			Q -= correction.getEntry(2);
		}
		System.out.println("E: " + E + ", L: " + L + ", Q: " + Q + ", Error = " + qDot().dotProduct(qDot()));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		InitialConditions ic = new InitialConditions(3.0, 3.0, PI / 2.0, -1.0);
		ic.solve();
		new KerrMotion(1.0, ic.a, 1.0, ic.E, ic.L, ic.Q, sqrt(ic.rMin * ic.rMax), PI / 2.0, 10.0, 0.001, 4).simulate();
	}

}