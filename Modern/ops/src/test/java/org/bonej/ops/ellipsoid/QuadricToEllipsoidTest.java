
package org.bonej.ops.ellipsoid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.Serializable;
import java.util.List;
import java.util.Optional;

import net.imagej.ImageJ;
import net.imagej.ops.special.function.BinaryFunctionOp;
import net.imagej.ops.special.function.Functions;
import net.imagej.ops.special.function.UnaryFunctionOp;
import net.imagej.ops.special.hybrid.BinaryHybridCFI1;
import net.imagej.ops.special.hybrid.Hybrids;

import org.bonej.ops.RotateAboutAxis;
import org.bonej.ops.SolveQuadricEq;
import org.junit.AfterClass;
import org.junit.Test;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Matrix4d;
import org.scijava.vecmath.Vector3d;

/**
 * Tests for {@link QuadricToEllipsoid}.
 *
 * @author Richard Domander
 */
public class QuadricToEllipsoidTest {

	private static final ImageJ IMAGE_J = new ImageJ();
	@SuppressWarnings("unchecked")
	private static UnaryFunctionOp<Matrix4d, Optional<Ellipsoid>> quadricToEllipsoid =
		(UnaryFunctionOp) Functions.unary(IMAGE_J.op(), QuadricToEllipsoid.class,
			Optional.class, Matrix4d.class);
	@SuppressWarnings("unchecked")
	private static BinaryFunctionOp<double[], Long, List<Vector3d>> ellipsoidPoints =
		(BinaryFunctionOp) Functions.binary(IMAGE_J.op(), EllipsoidPoints.class,
			List.class, new double[] { 1, 2, 3 }, 0);

	// TODO test ellipsoid "band"

	@Test
	public void testConeIsNotEllipsoid() throws Exception {
		//@formatter:off
		final Optional<Ellipsoid> coneOptional = quadricToEllipsoid.calculate(
				new Matrix4d(new double[] {
						1, 0, 0, 0,
						0, -1, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, -1 }));
		//@formatter:on

		assertFalse(coneOptional.isPresent());
	}

	@Test
	public void testCylinderIsNotEllipsoid() throws Exception {
		// Cylinder also causes SingularMatrixException which should be caught
		//@formatter:off
		final Optional<Ellipsoid> cylinderOptional = quadricToEllipsoid.calculate(
				new Matrix4d(new double[] {
						1, 0, 0, 0,
						0, 0, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, -1 }));
		//@formatter:on

		assertFalse(cylinderOptional.isPresent());
	}

	/**
	 * Tests fitting ellipsoid on a point cloud that forms a "band" around an
	 * ellipsoid surface, i.e. the radii of the points on the surface is scaled
	 * randomly
	 */
	@Test
	public void testEllipsoidBand() throws Exception {
		// SETUP
		final double[] radii = { 1, 2, 3 };

		// EXECUTE
		final List<Vector3d> points = ellipsoidPoints.calculate(radii, 1_000L);
		// The points are isotropically distributed on the ellipsoid surface, but
		// after the scaling they are not evenly distributed in space.
		points.forEach(p -> {
			final double scale = (2 * Math.random() - 1) * 0.05 + 1.0;
			p.scale(scale);
		});
		final Matrix4d quadric = (Matrix4d) IMAGE_J.op().run(SolveQuadricEq.class,
			points);
		final Optional<Ellipsoid> ellipsoidOptional = quadricToEllipsoid.calculate(
			quadric);

		// VERIFY
		assertTrue("Failed to fit transformed ellipsoid", ellipsoidOptional
			.isPresent());
		final Ellipsoid ellipsoid = ellipsoidOptional.get();
		assertTrue("Ellipsoid centre point is not within tolerance", new Vector3d(0,
			0, 0).epsilonEquals(ellipsoid.getCentroid(), 0.05));
		assertEquals(radii[0], ellipsoid.getA(), 0.025);
		assertEquals(radii[1], ellipsoid.getB(), 0.025);
		assertEquals(radii[2], ellipsoid.getC(), 0.025);
	}

	@Test
	public void testTransformedEllipsoid() throws Exception {
		// SETUP
		final Vector3d centroid = new Vector3d(1, 1, 1);
		final AxisAngle4d rotation = new AxisAngle4d(0, 0, 1, Math.PI / 4.0);
		final double[] radii = { 1, 2, 3 };

		// EXECUTE
		final List<Vector3d> points = ellipsoidPoints.calculate(radii, 1_000L);
		final BinaryHybridCFI1<Serializable, AxisAngle4d, Vector3d> rotate = Hybrids
			.binaryCFI1(IMAGE_J.op(), RotateAboutAxis.class, Vector3d.class,
				Vector3d.class, rotation);
		points.forEach(rotate::mutate);
		points.forEach(p -> p.add(centroid));
		final Matrix4d quadric = (Matrix4d) IMAGE_J.op().run(SolveQuadricEq.class,
			points);
		final Optional<Ellipsoid> ellipsoidOptional = quadricToEllipsoid.calculate(
			quadric);

		// VERIFY
		assertTrue("Failed to fit transformed ellipsoid", ellipsoidOptional
			.isPresent());
		final Ellipsoid transformedEllipsoid = ellipsoidOptional.get();
		assertTrue(transformedEllipsoid.getCentroid().epsilonEquals(centroid,
			1e-12));
		assertEquals(radii[0], transformedEllipsoid.getA(), 1e-12);
		assertEquals(radii[1], transformedEllipsoid.getB(), 1e-12);
		assertEquals(radii[2], transformedEllipsoid.getC(), 1e-12);
	}

	@Test
	public void testUnitSphere() throws Exception {
		//@formatter:off
		final Optional<Ellipsoid> sphereOptional = quadricToEllipsoid.calculate(
				new Matrix4d(new double[] {
						1, 0, 0, 0,
						0, 1, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, -1 }));
		//@formatter:on

		assertTrue("Failed to fit unit sphere", sphereOptional.isPresent());
		final Ellipsoid unitSphere = sphereOptional.get();
		assertEquals(1.0, unitSphere.getA(), 1e-12);
		assertEquals(1.0, unitSphere.getB(), 1e-12);
		assertEquals(1.0, unitSphere.getC(), 1e-12);
		assertTrue(unitSphere.getCentroid().epsilonEquals(new Vector3d(0, 0, 0),
			1e-12));
	}

	@AfterClass
	public static void oneTimeTearDown() {
		IMAGE_J.context().dispose();
	}
}
