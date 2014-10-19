/**
 * 
 */
package uk.me.doitto;

/*
 * #%L
 * blackhole4dclient
 * %%
 * Copyright (C) 2013 The Pond
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the 
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public 
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.json.simple.JSONObject;
import org.json.simple.JSONValue;

/**
 * @author ian
 * Particle trajectories in the Kerr spacetime and Boyer-Lindquist coordinates
 * <p>
 * This module loads initial conditions from JSON parameter files
 */
public class KerrSpaceTimeClient {
	
	double E, Lz, C;
	
	private void equatorial (double r, double a) {  // L and E for a circular orbit of r
		double sqrtR = Math.sqrt(r);
		double tmp = Math.sqrt(r * r - 3.0 * r + 2.0 * a * sqrtR);
		Lz = (r * r - 2.0 * a * sqrtR + a * a) / (sqrtR * tmp);
		E = (r * r - 2.0 * r + a * sqrtR) / (r * tmp);
	}

	private void constR (double M, double r, double a) {
		double B = 1 / (a * (M - r));
		double a2 = a * a;
		double r2 = r * r;
		double r4 = r2 * r2;
		double M2r = 2 * M * r;
		double E2 = E * E;
		double M2 = M * M;
		Lz = B * (M * E * (a2 - r2) + (a2 + r2 - M2r) * Math.sqrt(r2 * E2 + r * (M - r)));
		C = r2 * B * B * ((1 - E2) * r4 + (4 * M * E2 - 5 * M) * r2 * r + (8 - 5 * E2) * M2 * r2 +
				(2 * a2 * E2 - a2 - 4 * M2) * M * r + a2 * M2 + 2 * M * E * (a2 + r2 - M2r) * Math.sqrt(r * (M - r + r * E2)));
	}
	
	/**
	 * Read initial conditions from a JSON-formatted parameter file using Google's SimpleJSON library
	 * @param fileName the path to the file
	 * @return a KerrMotion instance
	 */
	public KerrMotion icJson (String fileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(fileName)));
		String data = "";
		String line = bufferedReader.readLine();
		while (line != null) {
			data += line;
			line = bufferedReader.readLine();
		}
		bufferedReader.close();
		JSONObject ic = (JSONObject)JSONValue.parse(data);
		Double M = (Double)ic.get("M");
		Double a = (Double)ic.get("a");
		Double mu = (Double)ic.get("mu");
		Double E = (Double)ic.get("E");
		Double Lz = (Double)ic.get("Lz");
		Double C = (Double)ic.get("C");
		Double K = (Double)ic.get("K");
		Double r = (Double)ic.get("r");
		Double theta = (Double)ic.get("theta");
		Double phi = (Double)ic.get("phi");
		Double time = (Double)ic.get("time");
		Double step = (Double)ic.get("step");
		int order = ((Long)ic.get("integratorOrder")).intValue();
		if (K != null) {  // circular polar if Lz == 0.0
			return new KerrMotion (M, a, mu, E, Lz, K - a * a * E * E, r, theta, phi, time, step, order);
		} else if (E == null && Lz == null && C == null) {  // constant radius
			equatorial(r, a);
			System.err.println("E: " + this.E + ", Lz: " + this.Lz);
			constR(M, r, a);
			System.err.println("E: " + this.E + ", Lz: " + this.Lz + ", CC: " + this.C);
			return new KerrMotion (M, a, mu, this.E, this.Lz, this.C, r, theta, phi, time, step, order);
		} else if (E == null && Lz == null && C == 0.0) {  // circular equatorial if C == 0.0
			equatorial(r, a);
			System.err.println("E: " + this.E + ", Lz: " + this.Lz);
			return new KerrMotion (M, a, mu, this.E, this.Lz, C, r, theta, phi, time, step, order);
		} else {
			return new KerrMotion (M, a, mu, E, Lz, C, r, theta, phi, time, step, order);
		}
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main (String[] args) throws IOException {
		if (args.length == 1) {
			new KerrSpaceTimeClient().icJson(args[0]).simulate();
		} else {
			System.err.println("Missing file name, giving up!");
		}
	}
}

