package uk.me.doitto;

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
    public void circleStable () {
		new KerrMotion(1.0, 0.0, 1.0, 0.9622504486493773, 4.0, 1.0, 12.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void circleUnStable () {
		new KerrMotion(1.0, 0.0, 1.0, 1.0, 4.0, 1.0, 4.0, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
    }

    @Test
    public void circleStablePrograde () {
		new KerrMotion(1.0, 1.0, 1.0, 0.9656115180853405, 4.0, 1.0, 14.109927135299376, Math.PI / 2.0, 0.0, 5000.0, 1.0, 4).simulate();
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