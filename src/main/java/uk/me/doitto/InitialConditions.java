/**
 * 
 */
package uk.me.doitto;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealVector;

/**
 * @author ian
 *
 */
public class InitialConditions {
	
	private double M = 1.0, mu = 1.0, E = 1.0, L = 2.0, Q = 0.0, a, rMin, rMax, thetaMax;
	
	/**
	 * 
	 */
	public InitialConditions(double rMin, double rMax, double thetaMax, double a) {
		this.rMin = rMin;
		this.rMax = rMax;
		this.thetaMax = thetaMax;
		this.a = a;
	}

	private double rDot (double r) {
		return ((r * r + a * a) * E - a * L) * ((r * r + a * a) * E - a * L) - (r * r - 2.0 * M * r + a * a) * (mu * mu * r * r + (L - a * E) + Q);
	}
	
	private double thetaDot (double theta) {
		return Q - cos(theta) * cos(theta) * (a * a * (mu * mu - E * E) + L * L / (sin(theta) * sin(theta)));
	}
	
	private RealVector qDot () {
		return new ArrayRealVector(new double[] { rDot(rMin), rDot(rMax), thetaDot(thetaMax) }, false);
	}
	
	private Array2DRowRealMatrix getJacobian () {
		return new Array2DRowRealMatrix(new double[][] {
			{ 2.0 * (rMin * rMin + a * a) * (E * (rMin * rMin + a * a) - a * L) + 2.0 * a * (L - a * E) * (rMin * rMin - 2.0 * M * rMin * rMin + a * a),
				-2.0 * a * (E * (rMin * rMin + a * a) - a * L) - 2.0 * (L - a * E) * (rMin * rMin - 2.0 * M * rMin + a * a),
				- rMin * rMin + 2.0 * M * rMin - a * a },
			{ 2.0 * (rMax * rMax + a * a) * (E * (rMax * rMax + a * a) - a * L) + 2.0 * a * (L - a * E) * (rMax * rMax - 2.0 * M * rMax * rMax + a * a),
				-2.0 * a * (E * (rMax * rMax + a * a) - a * L) - 2.0 * (L - a * E) * (rMax * rMax - 2.0 * M * rMax + a * a),
				- rMax * rMax + 2.0 * M * rMax - a * a }, 
			{ 2.0 * cos(thetaMax) * cos(thetaMax) * a * a * E, - (2.0 * cos(thetaMax) * cos(thetaMax) * L) / sin(thetaMax) * sin(thetaMax), 1.0 } }, false);
	}
	
	private void solve () {
		RealVector solution = new ArrayRealVector(new double[] { E, L, Q }, false);
		while (qDot().dotProduct(qDot()) > 1.0e-9) {
			RealVector newSolution = new LUDecomposition(getJacobian()).getSolver().solve(qDot());
//			RealVector newSolution = solution.subtract(solver.solve(qDot()));
//			solution = newSolution;
			E -= newSolution.getEntry(0);
			L -= newSolution.getEntry(1);
			Q -= newSolution.getEntry(2);
			System.out.println(newSolution);
		}
		System.out.println("E: " + E + ", L:" + L + ", Q: " + Q + ", Error = " + qDot().dotProduct(qDot()));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		RealMatrix jacobian = new Array2DRowRealMatrix(new double[][] {{ 2, 3, -2 }, { -1, 7, 6 }, { 4, -3, -5 } }, false);
//		DecompositionSolver solver = new LUDecomposition(jacobian).getSolver();
//
//		RealVector error = new ArrayRealVector(new double[] { 1, -2, 1 }, false);
//		RealVector solution = solver.solve(error);
//
//		System.out.println(solution);
		
		new InitialConditions(3.0, 6.0, PI / 4.0, 1.0).solve();
	}

}
