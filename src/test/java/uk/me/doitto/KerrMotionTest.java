package uk.me.doitto;

import org.junit.Ignore;
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
		new KerrMotion(1.0, 1.0, 1.0, 0.96, 1.98, 6.0, 12.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

	@Test
    public void polar () {
		new KerrMotion(1.0, 1.0, 1.0, 0.9558, 0.035991, 14.119546, 10.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void complex () {
		new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 2.4, 3.0, 12.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 6).simulate();
    }

    @Ignore
    @Test
    public void plummet () {
		new KerrMotion(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 100.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void circleStable () {
    	double r = 12.0;
    	double a = 0.0;
    	circular(r, a);
		new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void circleUnStable () {
    	double r = 4.0;
    	double a = 0.0;
    	circular(r, a);
		new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void circleStablePrograde () {
    	double r = 12.0;
    	double a = 1.0;
    	circular(r, a);
		new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void circleStableRetrograde () {
    	double r = 12.0;
    	double a = -1.0;
    	circular(r, a);
		new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void circleUnStableRetrograde () {
    	double r = 6.0;
    	double a = -1.0;
    	circular(r, a);
		new KerrMotion(1.0, a, 1.0, E, L, 0.0, r, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Test
    public void Polar1 () {
    	double K = 14.783;
    	double a = 0.8;
    	double E = 0.956;
		new KerrMotion(1.0, a, 1.0, E, 0.0, K - a * a * E * E, 10.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Ignore
    @Test
    public void first () {
		new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Ignore
    @Test
    public void second () {
		new KerrMotion(1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Ignore
    @Test
    public void third () {
		new KerrMotion(1.0, -1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Ignore
    @Test
    public void fourth () {
		new KerrMotion(1.0, -1.0, 1.0, 0.989352727272727, -4.683, 0.0, 12.201, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Ignore
    @Test
    public void fifth () {
		new KerrMotion(1.0, -1.0, 1.0, 1.0, 4.0, 0.0, 4.0, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }

    @Ignore
    @Test
    public void sixth () {
		new KerrMotion(1.0, -1.0, 1.0, 0.966, 4.066, 2.0, 17.488, Math.PI / 2.0, 0.0, 100.0, 0.01, 2).simulate();
    }
}