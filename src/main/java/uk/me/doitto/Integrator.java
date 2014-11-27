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

import static java.lang.Math.pow;

/**
 * @author ian
 * <p>
 * Symplectic integrators of various orders and types
 */
public enum Integrator {
	/**
	 * Stormer-Verlet 2nd-order
	 */
	SV2 (new double[] { 1.0 }),
	
	/**
	 * Yoshida 4th-order
	 */
	SV4 (new double[] { Constants.y, - Constants.y * Constants.CUBE_ROOT_2 }),
	
	/**
	 * Yoshida 6th-order
	 */
	SV6 (new double[] { 0.78451361047755726381949763, 0.23557321335935813368479318, -1.17767998417887100694641568, 1.31518632068391121888424973 }),
									
	/**
	 * Yoshida 8th-order
	 */
	SV8 (new double[] { 0.74167036435061295344822780, -0.40910082580003159399730010, 0.19075471029623837995387626, -0.57386247111608226665638773, 0.29906418130365592384446354, 0.33462491824529818378495798, 0.31529309239676659663205666, -0.79688793935291635401978884 }),
									
	/**
	 * Yoshida 10th-order
	 */
	SV10 (new double[] { 0.09040619368607278492161150, 0.53591815953030120213784983, 0.35123257547493978187517736, -0.31116802097815835426086544, -0.52556314194263510431065549, 0.14447909410225247647345695, 0.02983588609748235818064083, 0.17786179923739805133592238, 0.09826906939341637652532377, 0.46179986210411860873242126, -0.33377845599881851314531820, 0.07095684836524793621031152, 0.23666960070126868771909819, -0.49725977950660985445028388, -0.30399616617237257346546356, 0.05246957188100069574521612, 0.44373380805019087955111365 });

	private static class Constants {
		public final static double CUBE_ROOT_2 = pow(2.0, 1.0 / 3.0);		
		public final static double y = 1.0 / (2.0 - CUBE_ROOT_2);
	}
	
	protected double[] gammas;
	
	Integrator (double[] coefficients) {
		this.gammas = coefficients;
	}
	
	/**
	 * Basic 2nd-order Stormer-Verlet step which is composed into higher order methods
	 * @param bh the KerrMotion instance
	 * @param y composition coefficient
	 */
	protected final void sv (KerrMotion bh, double y) {
		double halfY = 0.5 * y;
		bh.updateQ(halfY);
		bh.updateP(y);
		bh.updateQ(halfY);
	}
	
	/**
	 * Perform one iteration step for the configured integrator
	 * @param bh the KerrMotion object reference, for passing through to the Q & P update methods
	 */
	void solve (KerrMotion bh) {
		int tmp = gammas.length - 1;
		for (int i = 0; i < tmp; i++) {
			sv(bh, gammas[i]);
		}
		for (int i = tmp; i >= 0; i--) {
			sv(bh, gammas[i]);
		}
	}
}
