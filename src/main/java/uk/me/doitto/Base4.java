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

public class Base4 extends Integrator {

	public Base4(IModel model, double[] compWeight) {
		super(model, compWeight);
        double cbrt2 = pow(2.0, (1.0 / 3.0));
        this.baseCoeff = new double[]{ 0.5 / (2.0 - cbrt2), 1.0 / (2.0 - cbrt2), 0.5 * (1.0 - cbrt2) / (2.0 - cbrt2), - cbrt2 / (2.0 - cbrt2) };
	}

	@Override
	protected void integrate(double w) {
	    model.pUp(w * baseCoeff[0]);
	    model.qUp(w * baseCoeff[1]);
	    model.pUp(w * baseCoeff[2]);
	    model.qUp(w * baseCoeff[3]);
	    model.pUp(w * baseCoeff[2]);
	    model.qUp(w * baseCoeff[1]);
	    model.pUp(w * baseCoeff[0]);
	}
}
