package uk.me.doitto;

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
		double sqrtR = Math.sqrt(r);
		double tmp = Math.sqrt(r * r - 3.0 * r + 2.0 * a * sqrtR);
		L = (r * r - 2.0 * a * sqrtR + a * a) / (sqrtR * tmp);
		E = (r * r - 2.0 * r + a * sqrtR) / (r * tmp);
	}

	@Test
    public void polar1 () {
		assertTrue(new KerrMotion(1.0, 1.0, 1.0, 0.9558, 0.035991, 14.119546, 10.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0);
	}

    @Test
    public void polar2 () {
    	double K = 14.783;
    	double a = 0.8;
    	double E = 0.956;
    	assertTrue(new KerrMotion(1.0, a, 1.0, E, 0.0, K - a * a * E * E, 10.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

   @Test
    public void complex () {
    	assertTrue(new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 2.4, 3.0, 12.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 6).simulate() < 1.0e-3);
    }

   @Test
   public void interesting () {
   		assertTrue(new KerrMotion(1.0, 1.0, 1.0, 0.96, 1.98, 6.8, 12.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 4).simulate() < 1.0e-3);
   }

    @Test
    public void plummetAZero () {
       	assertTrue(new KerrMotion(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 20.0, Math.PI / 2.0, 0.0, 100.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void plummetAPlus () {
    	assertTrue(new KerrMotion(1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 20.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void plummetAMinus () {
    	assertTrue(new KerrMotion(1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 20.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void circleStable () {
    	double r = 12.0;
    	double a = 0.0;
    	circular(r, a);
    		assertTrue(new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 10.0, 0.001, 2).simulate() < 1.0e-3);
    }

    @Test
    public void circleUnStable () {
    	double r = 4.0;
    	double a = 0.0;
    	circular(r, a);
    	assertTrue(new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 10.0, 0.001, 2).simulate() < 1.0e-3);
    }

    @Test
    public void circleStablePrograde () {
    	double r = 12.0;
    	double a = 1.0;
    	circular(r, a);
    	assertTrue(new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 10.0, 0.001, 2).simulate() < 1.0e-3);
    }

    @Test
    public void circleStableRetrograde () {
    	double r = 12.0;
    	double a = -1.0;
    	circular(r, a);
    	assertTrue(new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 10.0, 0.001, 2).simulate() < 1.0e-3);
    }

    @Test
    public void circleUnStableRetrograde () {
    	double r = 6.0;
    	double a = -1.0;
    	circular(r, a);
    	assertTrue(new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 10.0, 0.001, 2).simulate() < 1.0e-3);
    }

    @Test
    public void retrograde () {
    	assertTrue(new KerrMotion(1.0, -1.0, 1.0, 1.0, 4.0, 0.0, 4.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void messyInfall () {
       	assertTrue(new KerrMotion(1.0, 0.0, 1.0, 0.95, 1.5, 7.0, 14.64473, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void nearInfallAPlus () {
    	assertTrue(new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 2.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void infallAMinus () {
    	assertTrue(new KerrMotion(1.0, -1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 10.0, 0.001, 6).simulate() < 1.0e-3);
    }

    @Test
    public void precessingAMinus () {
    	assertTrue(new KerrMotion(1.0, -1.0, 1.0, 0.989352727272727, -4.683, 0.0, 12.201, Math.PI / 2.0, 0.0, 10.0, 0.001, 8).simulate() < 1.0e-3);
    }

    @Test
    public void spherical () {
    	assertTrue(new KerrMotion(1.0, 0.95, 1.0, 0.9, 1.707058, 4.021303, 4.0, Math.PI / 2.0, 0.0, 10.0, 0.01, 6).simulate() < 1.0e-3);
    }
}