package uk.me.doitto;

import static java.lang.Math.sqrt;
import static java.lang.Math.PI;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Tests for {@link KerrMotion}.
 *
 * @author ian
 */
@RunWith(JUnit4.class)
public class KerrMotionTest {
	
	private double L, E;

	private void circular (double r, double a) {  // L and E for a circular orbit of r
		double sqrtR = sqrt(r);
		double tmp = sqrt(r * r - 3.0 * r + 2.0 * a * sqrtR);
		L = (r * r - 2.0 * a * sqrtR + a * a) / (sqrtR * tmp);
		E = (r * r - 2.0 * r + a * sqrtR) / (r * tmp);
	}

	/**
	 * light 1-6: Edward Teo, Spherical photon orbits around a Kerr black hole
	 */
	@Test
    public void light1 () {
		double error = new KerrMotion(1.0, 1.0, 0.0, 1.0, 0.0, 22.31370849898476, 2.414213562373095, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-8);
	}

	@Test
    public void light2 () {
		double error = new KerrMotion(1.0, 1.0, 0.0, 1.0, -1.0, 25.856406460551018, 2.732050807568877, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-8);
	}

	@Test
    public void light3 () {
		double error = new KerrMotion(1.0, 1.0, 0.0, 1.0, -2.0, 27.0, 3.0, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-8);
	}

	@Test
    public void light4 () {
		double error = new KerrMotion(1.0, 1.0, 0.0, 1.0, -6.0, 9.627416997969522, 3.8284271247461903, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-8);
	}

	@Test
    public void light5 () {
		double error = new KerrMotion(1.0, 1.0, 0.0, 1.0, 1.0, 16.0, 2.0, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-8);
	}

	@Test
    public void light6 () {
		double error = new KerrMotion(1.0, 1.0, 0.0, 1.0, 1.999, 3.2590, 1.0316, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-4);
	}

	@Test
    public void polar1 () {
		double error = new KerrMotion(1.0, 1.0, 1.0, 0.9558, 0.035991, 14.119546, 10.0, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 0.2);
	}

    @Test
    public void polar2 () {
    	double K = 14.783;
    	double a = 0.8;
    	double E = 0.956;
		double error = new KerrMotion(1.0, a, 1.0, E, 0.0, K - a * a * E * E, 10.0, PI / 2.0, 10.0, 0.001, 8).simulate();
		assertTrue(error < 1.0e-8);
    }

   @Test
    public void complex () {
	   double error = new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 2.4, 3.0, 12.0, PI / 2.0, 10.0, 0.001, 6).simulate();
	   assertTrue(error < 1.0e-6);
    }

   @Test
   public void interesting () {
	   double error = new KerrMotion(1.0, 1.0, 1.0, 0.96, 1.98, 6.8, 12.0, PI / 2.0, 10.0, 0.001, 4).simulate();
   		assertTrue(error < 1.0e-4);
   }

    @Test
    public void plummetAZero () {
    	double error = new KerrMotion(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 20.0, PI / 2.0, 100.0, 0.001, 8).simulate();
       	assertTrue(error < 1.0e-8);
    }

    @Test
    public void plummetAPlus () {
    	double error = new KerrMotion(1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 20.0, PI / 2.0, 10.0, 0.001, 8).simulate();
    	assertTrue(error < 1.0e-7);
    }

    @Test
    public void plummetAMinus () {
    	double error = new KerrMotion(1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 20.0, PI / 2.0, 10.0, 0.001, 8).simulate();
    	assertTrue(error < 1.0e-7);
    }

    @Test
    public void circleStable () {
    	double r = 12.0;
    	double a = 0.0;
    	circular(r, a);
    	double error = new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, PI / 2.0, 10.0, 0.001, 2).simulate();
    	assertTrue(error < 1.0e-8);
    }

    @Test
    public void circleUnStable () {
    	double r = 4.0;
    	double a = 0.0;
    	circular(r, a);
    	double error = new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, PI / 2.0, 10.0, 0.001, 2).simulate();
    	assertTrue(error < 1.0e-24);
    }

    @Test
    public void circleStablePrograde () {
    	double r = 12.0;
    	double a = 1.0;
    	circular(r, a);
    	double error = new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, PI / 2.0, 10.0, 0.001, 2).simulate();
    	assertTrue(error < 1.0e-7);
    }

    @Test
    public void circleStableRetrograde () {
    	double r = 12.0;
    	double a = -1.0;
    	circular(r, a);
    	double error = new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, PI / 2.0, 10.0, 0.001, 2).simulate();
    	assertTrue(error < 1.0e-7);
    }

    @Test
    public void circleUnStableRetrograde () {
    	double r = 6.0;
    	double a = -1.0;
    	circular(r, a);
    	double error = new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, PI / 2.0, 10.0, 0.001, 2).simulate();
    	assertTrue(error < 1.0e-21);
    }

    @Test
    public void retrograde () {
    	double error = new KerrMotion(1.0, -1.0, 1.0, 1.0, 4.0, 0.0, 4.0, PI / 2.0, 10.0, 0.001, 8).simulate();
    	assertTrue(error < 1.0e-10);
    }

    @Test
    public void messyInfall () {
    	double error = new KerrMotion(1.0, 0.0, 1.0, 0.95, 1.5, 7.0, 14.64473, PI / 2.0, 10.0, 0.001, 8).simulate();
       	assertTrue(error < 1.0e-8);
    }

    @Test
    public void nearInfallAPlus () {
    	double error = new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 2.0, 1.0, 12.0, PI / 2.0, 10.0, 0.001, 8).simulate();
    	assertTrue(error < 1.0e-6);
    }

    @Test
    public void infallAMinus () {
    	double error = new KerrMotion(1.0, -1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 12.0, PI / 2.0, 10.0, 0.001, 6).simulate();
    	assertTrue(error < 1.0e-8);
    }

    @Test
    public void precessingAMinus () {
    	double error = new KerrMotion(1.0, -1.0, 1.0, 0.989352727272727, -4.683, 0.0, 12.201, PI / 2.0, 10.0, 0.001, 8).simulate();
    	assertTrue(error < 1.0e-4);
    }

    @Test
    public void spherical1 () {
    	double error = new KerrMotion(1.0, 0.95, 1.0, 0.9, 1.707058, 4.021303, 4.0, PI / 2.0, 10.0, 0.01, 6).simulate();
    	assertTrue(error < 1.0e-7);
    }

    @Test
    public void spherical2 () {
    	double error = new KerrMotion(1.0, 1.0, 1.0, 0.9601, -2.667245, 9.381218, 10.0, PI / 2.0, 10.0, 0.01, 6).simulate();
    	assertTrue(error < 1.0e-5);
    }
}