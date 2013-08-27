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

    @Test
    public void Polar () {
		new KerrMotion(1.0, 1.0, 1.0, 0.96, 1.98, 6.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 2).simulate();
    }

    @Test
    public void Complex () {
		new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 2.4, 3.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 6).simulate();
    }

    @Ignore
    @Test
    public void Plummet () {
		new KerrMotion(1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 100.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void circleStable () {
		new KerrMotion(1.0, 0.0, 1.0, 0.9622504486493773, 4.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 8).simulate();
    }

    @Test
    public void circleUnStable () {
		new KerrMotion(1.0, 0.0, 1.0, 1.0, 4.0, 0.0, 4.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 10).simulate();
    }

    @Test
    public void circleStablePrograde () {
		new KerrMotion(1.0, 1.0, 1.0, 0.959723537290, 3.717928598, 0.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 2).simulate();
    }

    @Test
    public void circleStableRetrograde () {
		new KerrMotion(1.0, -1.0, 1.0, 0.965969668, 4.362473, 0.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 2).simulate();
    }

    @Test
    public void circleUnStableRetrograde () {
		new KerrMotion(1.0, -1.0, 1.0, 0.992324568078, 4.725799115, 0.0, 6.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 2).simulate();
    }

    @Test
    public void first () {
		new KerrMotion(1.0, 1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void second () {
		new KerrMotion(1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void third () {
		new KerrMotion(1.0, -1.0, 1.0, 0.962250448649377, 0.6 * 4.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void fourth () {
		new KerrMotion(1.0, -1.0, 1.0, 0.989352727272727, -4.683, 0.0, 12.201, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void fifth () {
		new KerrMotion(1.0, -1.0, 1.0, 1.0, 4.0, 0.0, 4.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void sixth () {
		new KerrMotion(1.0, -1.0, 1.0, 0.966, 4.066, 2.0, 17.488, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }
}