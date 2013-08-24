package uk.me.doitto;

import static org.junit.Assert.fail;
import static uk.me.doitto.Integrator.STORMER_VERLET_10;
import static uk.me.doitto.Integrator.STORMER_VERLET_2;
import static uk.me.doitto.Integrator.STORMER_VERLET_4;
import static uk.me.doitto.Integrator.STORMER_VERLET_6;
import static uk.me.doitto.Integrator.STORMER_VERLET_8;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

public class IntegratorTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Ignore
	@Test
	public void testSympBase() {
	}

	@Test
	public void testInit() {
		STORMER_VERLET_2.init();
		STORMER_VERLET_4.init();
		STORMER_VERLET_6.init();
		STORMER_VERLET_8.init();
		STORMER_VERLET_10.init();
	}

	@Ignore
	@Test
	public void testSolve() {
		fail("Not yet implemented");
	}

}
