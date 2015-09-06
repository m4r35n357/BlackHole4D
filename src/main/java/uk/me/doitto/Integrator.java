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


/**
 * Integrator superclass with an ISymplectic factory method
 */
public abstract class Integrator implements ISymplectic {

    private double[] compWeight;
    private int wRange;
    protected IModel model;
    protected double[] baseCoeff;

    protected Integrator (IModel model, double[] compWeight) {
        this.compWeight = compWeight;
        this.wRange = compWeight.length;
        this.model = model;
    }

	public static ISymplectic getIntegrator (IModel model, String order) {
		ISymplectic integrator = null;
	    switch (order) {
	        case "sb2":  // second order, basic
	            integrator = new Base2(model, new double[]{ 1.0 });
	            break;
	        case "sc4":  // fourth order, composed
				double cbrt2 = pow(2.0, (1.0 / 3.0));
	            integrator = new Base2(model, new double[]{ 1.0 / (2.0 - cbrt2), - cbrt2 / (2.0 - cbrt2), 1.0 / (2.0 - cbrt2) });
	            break;
	        case "sb4":  // fourth order, basic
	            integrator = new Base4(model, new double[]{ 1.0 });
	            break;
	        case "sc6":  // sixth order, composed
	            double fthrt2 = pow(2.0, (1.0 / 5.0));
	            integrator = new Base4(model, new double[]{ 1.0 / (2.0 - fthrt2), - fthrt2 / (2.0 - fthrt2), 1.0 / (2.0 - fthrt2) });
	            break;
	    }
        return integrator;
	}

    protected abstract void integrate (double w);

    public void compose () {
		for (int i = 0; i < wRange; i++) {
            integrate(compWeight[i] * model.getH());
        }
    }
}

