/*
Copyright (c) 2014, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */
package uk.me.doitto;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

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
		boolean nonsingular = abs(rMax - rMin) > 2.0 * tolerance;
		this.r0 = nonsingular ? rMin: rMin - tolerance;
		this.r1 = nonsingular ? rMax: rMax + tolerance;
		this.th0 = thetaMin > 0.01 ? thetaMin:  0.01;
		switch (a) {
			case RETROGRADE: this.a = -1.0; break;
			case ZERO: this.a = 0.0; break;
			case PROGRADE: this.a = 1.0; break;
		default: this.a = 0.0; break;
		}
		this.factorL = factorL;
		switch (i) {
			case SV2: this.order = 2; break;
			case SV4: this.order = 4; break;
			case SV6: this.order = 6; break;
			case SV8: this.order = 8; break;
			case SV10: this.order = 10; break;
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
		while (qDot().dotProduct(qDot()) > 1.0e-21) {
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
		InitialConditions ic = new InitialConditions(Trajectory.PARTICLE, 5.9, 5.9, PI / 2, Spin.ZERO, 1.0, Integrator.SV8);
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
