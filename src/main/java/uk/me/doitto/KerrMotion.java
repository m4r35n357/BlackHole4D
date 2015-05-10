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

import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.log10;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static uk.me.doitto.Integrator.SV10;
import static uk.me.doitto.Integrator.SV2;
import static uk.me.doitto.Integrator.SV4;
import static uk.me.doitto.Integrator.SV6;
import static uk.me.doitto.Integrator.SV8;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.json.simple.JSONObject;
import org.json.simple.JSONValue;

/**
 * @author ian
 * Particle and photon trajectories in the Kerr spacetime and in Boyer-Lindquist coordinates
 */
public final class KerrMotion {
	
	private final double M, a, horizon, mu2, E, E2, L, L2, Q, T, ts, a2, aE, a2E, aL, l_ae2, a2mu2_E2, nf = 1.0e-18; // constants for this spacetime
	
	private double r2, ra2, sth, cth, sth2, cth2, delta, R, P1, P2, THETA, TH;  // intermediate variables
	
	private double mino, tau, t, r, th, ph, tDot, rDot, thDot, phDot, eCum, e, eR, eTh; // coordinates etc.
	
	private Integrator integrator = SV2;
	
	/**
	 * Constructor, constants and initial conditions
	 */
	public KerrMotion (double bhMass, double spin, double pMass, double energy, double zAngMom, double CC, double r0, double th0, double duration, double timestep, int order) {
		M = bhMass;
		a = spin;
		a2 = a * a;
		horizon = M * (1.0 + sqrt(1.0 - a2));
		mu2 = pMass * pMass;
		E = energy;
		E2 = E * E;
		aE = a * E;
		a2E = a2 * E;
		L = zAngMom;
		L2 = L * L;
		aL = a * L;
		l_ae2 = (L - a * E) * (L - a * E);
		a2mu2_E2 = a2 * (mu2 - E2);
		Q = CC;
		r = r0;
		th = th0;
		T = duration;
		ts = timestep;
		switch (order) {
			case 2: integrator = SV2; break;
			case 4: integrator = SV4; break;
			case 6: integrator = SV6; break;
			case 8: integrator = SV8; break;
			case 10: integrator = SV10; break;
		}
	}

	private double clamp (double potential) {
		return potential >= 0.0 ? potential : 0.0;
	}
	
	private void updateIntermediates () {
		r2 = r * r;
		ra2 = r2 + a2;
		sth = sin(th);
		cth = cos(th);
		sth2 = sth * sth;
		cth2 = cth * cth;
		delta = ra2 - 2.0 * M * r;
		P1 = ra2 * E - aL;  // MTW eq.33.33b, ignoring charge term
		P2 = mu2 * r2 + l_ae2 + Q;
		R = P1 * P1 - delta * P2;  // MTW eq.33.33c
		TH = a2mu2_E2 + L2 / sth2;
		THETA = Q - cth2 * TH;  // MTW eq.33.33a
	}
	
	private void errors () {
		double e_r = abs(rDot * rDot - clamp(R)) / 2.0;
		double e_th = abs(thDot * thDot - clamp(THETA)) / 2.0;
		eR = 10.0 * log10(e_r >= nf ? e_r : nf);
		eTh = 10.0 * log10(e_th >= nf ? e_th : nf);
		e =  10.0 * log10(e_r + e_th >= nf ? e_r + e_th: nf);
		eCum += e_r + e_th;
	}
	
	private void update_t_phi_Dot () {
		tDot = ra2 * P1 / delta + aL - a2E * sth2;  // MTW eq.33.32d
		phDot = a * P1 / delta - aE + L / sth2;  // MTW eq.33.32c
	}
	
	private void update_t_phi () {
		update_t_phi_Dot();
		t += ts * tDot;
		ph += ts * phDot;
	}
	
	void updateQ (double c) {  // dH/dX
		r += c * ts * rDot;
		th += c * ts * thDot;
		updateIntermediates();
	}
	
	void updateP (double c) {  // dH/dXdot
		rDot += c * ts * (2.0 * r * E * P1 - P2 * (r - M) - mu2 * r * delta);  // dR/dr see Maxima file maths.wxm, "My Equations (Mino Time)"
		thDot += c * ts * (cth * sth * TH + L2 * cth2 * cth / (sth2 * sth));  // dTheta/dtheta see Maxima file maths.wxm, "My Equations (Mino Time)"
	}
	
	public double simulate () {
		updateIntermediates();
		rDot = - sqrt(clamp(R));  // MTW eq.33.32b
		thDot = - sqrt(clamp(THETA));  // MTW eq.33.32a
		update_t_phi_Dot();
		do {
			double ra = sqrt(ra2);
			double sigma = (r2 + a2 * cth2);
			errors();
			System.out.printf("{\"mino\":%.9e, \"tau\":%.9e, \"E\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, \"EC\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tDot\":%.9e, \"rDot\":%.9e, \"thDot\":%.9e, \"phDot\":%.9e, \"x\":%.9e, \"y\":%.9e, \"z\":%.9e}%n",
					mino, tau, e, eR, eTh, 10.0 * log10(eCum >= nf ? eCum : nf), t, r, th, ph, tDot / sigma, rDot / sigma, thDot / sigma, phDot / sigma, ra * sth * cos(ph), ra * sth * sin(ph), r * cth);
			update_t_phi();  // Euler
			integrator.solve(this);
			mino += ts;
			tau += ts * sigma;
		} while (r > horizon && mino <= T);  // outside horizon and in proper time range
		return eCum;
	}
	
	/**
	 * Read initial conditions from a JSON-formatted parameter file using Google's SimpleJSON library
	 * 
	 * @param args[0] the path to the file, if no args then read from stdin
	 * @throws IOException 
	 */
	public static void main (String[] args) throws IOException {
		BufferedReader bufferedReader;
		if (args.length == 1) {
			bufferedReader = new BufferedReader(new FileReader(new File(args[0])));
		} else {
			bufferedReader = new BufferedReader(new InputStreamReader(System.in));
		}
		StringBuilder data = new StringBuilder();
		String line = bufferedReader.readLine();
		while (line != null) {
			data.append(line);
			line = bufferedReader.readLine();
		}
		bufferedReader.close();
		JSONObject ic = (JSONObject)JSONValue.parse(data.toString());
		new KerrMotion ((double)ic.get("M"), (double)ic.get("a"), (double)ic.get("mu"), (double)ic.get("E"), (double)ic.get("Lz"), (double)ic.get("C"),
			(double)ic.get("r"), (double)ic.get("theta"), (double)ic.get("time"), (double)ic.get("step"), ((Long)ic.get("integratorOrder")).intValue()).simulate();
	}
}
