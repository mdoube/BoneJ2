package org.bonej.plugins;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.concurrent.atomic.AtomicInteger;

import org.bonej.geometry.FitEllipsoid;
import org.bonej.util.Multithreader;
import org.scijava.vecmath.Point3f;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import customnode.CustomTriangleMesh;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import marchingcubes.MCTriangulator;

public class ParticleAnalyser {
	
	/** Foreground value */
	static final int FORE = -1;
	/** Background value */
	static final int BACK = 0;
	 
	private static String sPhase = "";
	
	/** number of particle labels */
	private static int nParticles;
	
	/**
	 * Get the sizes of all the particles as a voxel count
	 *
	 * @param particleLabels particles in the image.
	 * @return particleSizes sizes of the particles.
	 */
	static long[] getParticleSizes(final int[][] particleLabels) {
		IJ.showStatus("Getting " + sPhase + " particle sizes");
		final int d = particleLabels.length;
		final int wh = particleLabels[0].length;

		final int maxParticle = nParticles;
		
		//make a list of all the particle sizes with 
		//index = particle value
		//need to handle the ID offset for the chunks
		AtomicInteger an = new AtomicInteger(0);
		final long[][] partSizes = new long[d][];
		
		final Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(() -> {
				for (int z = an.getAndIncrement(); z < d; z = an.getAndIncrement()) {
					final long[] particleSizes = new long[maxParticle + 1];
					final int[] slice = particleLabels[z];
					for (int i = 0; i < wh; i++) {
						final int label = slice[i];
						//hack to avoid AIOOB
						if (label <= maxParticle)
							particleSizes[label]++;
					}
					partSizes[z] = particleSizes;
				}
			});
		}
		Multithreader.startAndJoin(threads);
		
		final long[] particleSizes = new long[maxParticle + 1];
		for (int i = 0; i <= maxParticle; i++) {
			long partSum = 0;
			for (int z = 0; z < d; z++)
				partSum += partSizes[z][i];
			particleSizes[i] = partSum;			
		}
		
		IJ.showStatus("Finished calculating particle sizes");
		return particleSizes;
	}
	
	/**
	 * Calculate surface area of the isosurface
	 *
	 * @param points in 3D triangle mesh
	 * @return surface area
	 */
	private static double getSurfaceArea(final List<Point3f> points) {
		if (points == null) {
			return 0;
		}
		double sumArea = 0;
		final int nPoints = points.size();
		final Point3f origin = new Point3f(0.0f, 0.0f, 0.0f);
		for (int n = 0; n < nPoints; n += 3) {
			IJ.showStatus("Calculating surface area...");
			final Point3f cp = crossProduct(points.get(n), points.get(n + 1), points
				.get(n + 2));
			final double deltaArea = 0.5 * cp.distance(origin);
			sumArea += deltaArea;
		}
		return sumArea;
	}

	static double[] getSurfaceAreas(
		final Collection<List<Point3f>> surfacePoints)
	{
		return surfacePoints.stream().mapToDouble(ParticleAnalyser::getSurfaceArea)
			.toArray();
	}

	static ArrayList<List<Point3f>> getSurfacePoints(final ImagePlus imp,
		final int[][] particleLabels, final int[][] limits, final int resampling)
	{
		final Calibration cal = imp.getCalibration();
		final ArrayList<List<Point3f>> surfacePoints = new ArrayList<>();
		final boolean[] channels = { true, false, false };
		for (int p = 0; p < nParticles; p++) {
			IJ.showStatus("Getting surface meshes...");
			IJ.showProgress(p, nParticles);
			if (p > 0) {
				final ImagePlus binaryImp = getBinaryParticle(p, imp, particleLabels,
					limits, resampling);
				// noinspection TypeMayBeWeakened
				final MCTriangulator mct = new MCTriangulator();
				@SuppressWarnings("unchecked")
				final List<Point3f> points = mct.getTriangles(binaryImp, 128, channels,
					resampling);

				final double xOffset = (limits[p][0] - 1) * cal.pixelWidth;
				final double yOffset = (limits[p][2] - 1) * cal.pixelHeight;
				final double zOffset = (limits[p][4] - 1) * cal.pixelDepth;
				for (final Point3f point : points) {
					point.x += xOffset;
					point.y += yOffset;
					point.z += zOffset;
				}
				surfacePoints.add(points);
				if (points.isEmpty()) {
					IJ.log("Particle " + p + " resulted in 0 surface points");
				}
			}
			else {
				surfacePoints.add(null);
			}
		}
		return surfacePoints;
	}

	static double[] getSurfaceVolume(
		final Collection<List<Point3f>> surfacePoints)
	{
		return surfacePoints.stream().mapToDouble(p -> {
			if (p == null) {
				return 0;
			}
			final CustomTriangleMesh mesh = new CustomTriangleMesh(p);
			return Math.abs(mesh.getVolume());
		}).toArray();
	}

	static double[] getVolumes(final ImagePlus imp,
		final long[] particleSizes)
	{
		final Calibration cal = imp.getCalibration();
		final double voxelVolume = cal.pixelWidth * cal.pixelHeight *
			cal.pixelDepth;
		final int nLabels = particleSizes.length;
		final double[] particleVolumes = new double[nLabels];
		for (int i = 0; i < nLabels; i++) {
			particleVolumes[i] = voxelVolume * particleSizes[i];
		}
		return particleVolumes;
	}

	static boolean isRadiiValid(final double[] radii) {
		for (int r = 0; r < 3; r++) {
			if (Double.isNaN(radii[r])) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Calculate the cross product of 3 Point3f's, which describe two vectors
	 * joined at the tails. Can be used to find the plane / surface normal of a
	 * triangle. Half of its magnitude is the area of the triangle.
	 *
	 * @param point0 both vectors' tails
	 * @param point1 vector 1's head
	 * @param point2 vector 2's head
	 * @return cross product vector
	 */
	private static Point3f crossProduct(final Point3f point0,
		final Point3f point1, final Point3f point2)
	{
		final double x1 = point1.x - point0.x;
		final double y1 = point1.y - point0.y;
		final double z1 = point1.z - point0.z;
		final double x2 = point2.x - point0.x;
		final double y2 = point2.y - point0.y;
		final double z2 = point2.z - point0.z;
		final Point3f crossVector = new Point3f();
		crossVector.x = (float) (y1 * z2 - z1 * y2);
		crossVector.y = (float) (z1 * x2 - x1 * z2);
		crossVector.z = (float) (x1 * y2 - y1 * x2);
		return crossVector;
	}

	/**
	 * create a binary ImagePlus containing a single particle and which 'just
	 * fits' the particle
	 *
	 * @param p The particle ID to get
	 * @param imp original image, used for calibration
	 * @param particleLabels work array of particle labels
	 * @param limits x,y and z limits of each particle
	 * @param padding amount of empty space to pad around each particle
	 * @return a cropped single particle image.
	 */
	private static ImagePlus getBinaryParticle(final int p, final ImagePlus imp,
		final int[][] particleLabels, final int[][] limits, final int padding)
	{

		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int xMin = Math.max(0, limits[p][0] - padding);
		final int xMax = Math.min(w - 1, limits[p][1] + padding);
		final int yMin = Math.max(0, limits[p][2] - padding);
		final int yMax = Math.min(h - 1, limits[p][3] + padding);
		final int zMin = Math.max(0, limits[p][4] - padding);
		final int zMax = Math.min(d - 1, limits[p][5] + padding);
		final int stackWidth = xMax - xMin + 1;
		final int stackHeight = yMax - yMin + 1;
		final int stackSize = stackWidth * stackHeight;
		final ImageStack stack = new ImageStack(stackWidth, stackHeight);
		for (int z = zMin; z <= zMax; z++) {
			final byte[] slice = new byte[stackSize];
			int i = 0;
			for (int y = yMin; y <= yMax; y++) {
				final int sourceIndex = y * w;
				for (int x = xMin; x <= xMax; x++) {
					if (particleLabels[z][sourceIndex + x] == p) {
						slice[i] = (byte) 0xFF;
					}
					i++;
				}
			}
			stack.addSlice(imp.getStack().getSliceLabel(z + 1), slice);
		}
		final ImagePlus binaryImp = new ImagePlus("Particle_" + p, stack);
		final Calibration cal = imp.getCalibration();
		binaryImp.setCalibration(cal);
		return binaryImp;
	}
	/**
	 * Get the centroids of all the particles in real units
	 *
	 * @param imp an image.
	 * @param particleLabels particles in the image.
	 * @param particleSizes sizes of the particles
	 * @return double[][] containing all the particles' centroids
	 */
	static double[][] getCentroids(final ImagePlus imp,
		final int[][] particleLabels, final long[] particleSizes)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final double[][] sums = new double[nParticles][3];
		for (int z = 0; z < d; z++) {
			for (int y = 0; y < h; y++) {
				final int index = y * w;
				for (int x = 0; x < w; x++) {
					final int particle = particleLabels[z][index + x];
					sums[particle][0] += x;
					sums[particle][1] += y;
					sums[particle][2] += z;
				}
			}
		}
		final Calibration cal = imp.getCalibration();
		final double[][] centroids = new double[nParticles][3];
		for (int p = 0; p < nParticles; p++) {
			centroids[p][0] = cal.pixelWidth * sums[p][0] / particleSizes[p];
			centroids[p][1] = cal.pixelHeight * sums[p][1] / particleSizes[p];
			centroids[p][2] = cal.pixelDepth * sums[p][2] / particleSizes[p];
		}
		return centroids;
	}

	static EigenvalueDecomposition[] getEigens(final ImagePlus imp,
		final int[][] particleLabels, final double[][] centroids)
	{
		final Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;
		final double voxVhVd = (vH * vH + vD * vD) / 12;
		final double voxVwVd = (vW * vW + vD * vD) / 12;
		final double voxVhVw = (vH * vH + vW * vW) / 12;
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final EigenvalueDecomposition[] eigens =
			new EigenvalueDecomposition[nParticles];
		final double[][] momentTensors = new double[nParticles][6];
		for (int z = 0; z < d; z++) {
			IJ.showStatus("Calculating particle moments...");
			IJ.showProgress(z, d);
			final double zVd = z * vD;
			for (int y = 0; y < h; y++) {
				final double yVh = y * vH;
				final int index = y * w;
				for (int x = 0; x < w; x++) {
					final int p = particleLabels[z][index + x];
					if (p > 0) {
						final double xVw = x * vW;
						final double dx = xVw - centroids[p][0];
						final double dy = yVh - centroids[p][1];
						final double dz = zVd - centroids[p][2];
						momentTensors[p][0] += dy * dy + dz * dz + voxVhVd; // Ixx
						momentTensors[p][1] += dx * dx + dz * dz + voxVwVd; // Iyy
						momentTensors[p][2] += dy * dy + dx * dx + voxVhVw; // Izz
						momentTensors[p][3] += dx * dy; // Ixy
						momentTensors[p][4] += dx * dz; // Ixz
						momentTensors[p][5] += dy * dz; // Iyz
					}
				}
			}
			for (int p = 1; p < nParticles; p++) {
				final double[][] inertiaTensor = new double[3][3];
				inertiaTensor[0][0] = momentTensors[p][0];
				inertiaTensor[1][1] = momentTensors[p][1];
				inertiaTensor[2][2] = momentTensors[p][2];
				inertiaTensor[0][1] = -momentTensors[p][3];
				inertiaTensor[0][2] = -momentTensors[p][4];
				inertiaTensor[1][0] = -momentTensors[p][3];
				inertiaTensor[1][2] = -momentTensors[p][5];
				inertiaTensor[2][0] = -momentTensors[p][4];
				inertiaTensor[2][1] = -momentTensors[p][5];
				final Matrix inertiaTensorMatrix = new Matrix(inertiaTensor);
				final EigenvalueDecomposition E = new EigenvalueDecomposition(
					inertiaTensorMatrix);
				eigens[p] = E;
			}
		}
		return eigens;
	}

	static Object[][] getEllipsoids(
		final Collection<List<Point3f>> surfacePoints)
	{
		final Object[][] ellipsoids = new Object[surfacePoints.size()][];
		int p = 0;
		for (final List<Point3f> points : surfacePoints) {
			if (points == null) {
				p++;
				continue;
			}
			final Iterator<Point3f> pointIter = points.iterator();
			final double[][] coOrdinates = new double[points.size()][3];
			int i = 0;
			while (pointIter.hasNext()) {
				final Point3f point = pointIter.next();
				coOrdinates[i][0] = point.x;
				coOrdinates[i][1] = point.y;
				coOrdinates[i][2] = point.z;
				i++;
			}
			try {
				ellipsoids[p] = FitEllipsoid.yuryPetrov(coOrdinates);
			}
			catch (final RuntimeException re) {
				IJ.log("Could not fit ellipsoid to surface " + p);
			}
			p++;
		}
		return ellipsoids;
	}

	/**
	 * Get the Euler characteristic of each particle
	 *
	 * @param imp an image.
	 * @param particleLabels particles of the image.
	 * @param limits limits of the particles.
	 * @return euler characteristic of each image.
	 */
	static double[][] getEulerCharacter(final ImagePlus imp,
		final int[][] particleLabels, final int[][] limits)
	{
		final Connectivity con = new Connectivity();
		final double[][] eulerCharacters = new double[nParticles][3];
		for (int p = 1; p < nParticles; p++) {
			final ImagePlus particleImp = getBinaryParticle(p, imp, particleLabels,
				limits, 1);
			final double euler = con.getSumEuler(particleImp);
			final double cavities = getNCavities(particleImp);
			// Calculate number of holes and cavities using
			// Euler = particles - holes + cavities
			// where particles = 1
			final double holes = cavities - euler + 1;
			final double[] bettis = { euler, holes, cavities };
			eulerCharacters[p] = bettis;
		}
		return eulerCharacters;
	}

	/**
	 * Get the Feret diameter of a surface. Uses an inefficient brute-force
	 * algorithm.
	 *
	 * @param particleSurfaces points of all the particles.
	 * @return Feret diameters of the surfaces.
	 */
	static double[] getFerets(
		final List<List<Point3f>> particleSurfaces)
	{
		final int nSurfaces = particleSurfaces.size();
		final double[] ferets = new double[nSurfaces];
		final ListIterator<List<Point3f>> it = particleSurfaces.listIterator();
		int i = 0;
		Point3f a;
		Point3f b;
		List<Point3f> surface;
		ListIterator<Point3f> ita;
		ListIterator<Point3f> itb;
		while (it.hasNext()) {
			IJ.showStatus("Finding Feret diameter...");
			IJ.showProgress(it.nextIndex(), nSurfaces);
			surface = it.next();
			if (surface == null) {
				ferets[i] = Double.NaN;
				i++;
				continue;
			}
			ita = surface.listIterator();
			while (ita.hasNext()) {
				a = ita.next();
				itb = surface.listIterator(ita.nextIndex());
				while (itb.hasNext()) {
					b = itb.next();
					ferets[i] = Math.max(ferets[i], a.distance(b));
				}
			}
			i++;
		}
		return ferets;
	}

	/**
	 * Get the mean and standard deviation of pixel values &gt;0 for each particle
	 * in a particle label work array
	 *
	 * @param imp Input image containing pixel values
	 * @param particleLabels workArray containing particle labels
	 * @param particleSizes array of particle sizes as pixel counts
	 * @return array containing mean, std dev and max pixel values for each
	 *         particle
	 */
	static double[][] getMeanStdDev(final ImagePlus imp,
		final int[][] particleLabels, final long[] particleSizes)
	{
		final int d = imp.getImageStackSize();
		final int wh = imp.getWidth() * imp.getHeight();
		final ImageStack stack = imp.getImageStack();
		final double[] sums = new double[nParticles];
		for (int z = 0; z < d; z++) {
			final float[] pixels = (float[]) stack.getPixels(z + 1);
			final int[] labelPixels = particleLabels[z];
			for (int i = 0; i < wh; i++) {
				final double value = pixels[i];
				if (value > 0) {
					sums[labelPixels[i]] += value;
				}
			}
		}
		final double[][] meanStdDev = new double[nParticles][3];
		for (int p = 1; p < nParticles; p++) {
			meanStdDev[p][0] = sums[p] / particleSizes[p];
		}

		final double[] sumSquares = new double[nParticles];
		for (int z = 0; z < d; z++) {
			final float[] pixels = (float[]) stack.getPixels(z + 1);
			final int[] labelPixels = particleLabels[z];
			for (int i = 0; i < wh; i++) {
				final double value = pixels[i];
				if (value > 0) {
					final int p = labelPixels[i];
					final double residual = value - meanStdDev[p][0];
					sumSquares[p] += residual * residual;
					meanStdDev[p][2] = Math.max(meanStdDev[p][2], value);
				}
			}
		}
		for (int p = 1; p < nParticles; p++) {
			meanStdDev[p][1] = Math.sqrt(sumSquares[p] / particleSizes[p]);
		}
		return meanStdDev;
	}

	private static int getNCavities(final ImagePlus imp) {
		final Object[] result = (new ParticleCounter()).getParticles(imp, BACK);
		final long[] particleSizes = (long[]) result[2];
		return particleSizes.length - 2;
	}

	/**
	 * Get the minimum and maximum x, y and z coordinates of each particle
	 *
	 * @param imp ImagePlus (used for stack size)
	 * @param particleLabels work array containing labelled particles
	 * @return int[][] containing x, y and z minima and maxima.
	 */
	static int[][] getParticleLimits(final ImagePlus imp,
		final int[][] particleLabels)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int[][] limits = new int[nParticles][6];
		for (int i = 0; i < nParticles; i++) {
			limits[i][0] = Integer.MAX_VALUE; // x min
			limits[i][1] = 0; // x max
			limits[i][2] = Integer.MAX_VALUE; // y min
			limits[i][3] = 0; // y max
			limits[i][4] = Integer.MAX_VALUE; // z min
			limits[i][5] = 0; // z max
		}
		for (int z = 0; z < d; z++) {
			for (int y = 0; y < h; y++) {
				final int index = y * w;
				for (int x = 0; x < w; x++) {
					final int i = particleLabels[z][index + x];
					limits[i][0] = Math.min(limits[i][0], x);
					limits[i][1] = Math.max(limits[i][1], x);
					limits[i][2] = Math.min(limits[i][2], y);
					limits[i][3] = Math.max(limits[i][3], y);
					limits[i][4] = Math.min(limits[i][4], z);
					limits[i][5] = Math.max(limits[i][5], z);
				}
			}
		}
		return limits;
	}
	
	/**
	 * Scans edge voxels and set all touching particles to background
	 *
	 * @param imp an image.
	 * @param particleLabels particles in the image.
	 */
	static void excludeOnEdges(final ImagePlus imp, final int[][] particleLabels,
		final byte[][] workArray)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final long[] particleSizes = getParticleSizes(particleLabels);
		final int nLabels = particleSizes.length;
		final int[] newLabel = new int[nLabels];
		for (int i = 0; i < nLabels; i++)
			newLabel[i] = i;

		// scan faces
		// top and bottom faces
		for (int y = 0; y < h; y++) {
			final int index = y * w;
			for (int x = 0; x < w; x++) {
				final int pt = particleLabels[0][index + x];
				if (pt > 0) newLabel[pt] = 0;
				final int pb = particleLabels[d - 1][index + x];
				if (pb > 0) newLabel[pb] = 0;
			}
		}

		// west and east faces
		for (int z = 0; z < d; z++) {
			for (int y = 0; y < h; y++) {
				final int pw = particleLabels[z][y * w];
				final int pe = particleLabels[z][y * w + w - 1];
				if (pw > 0) newLabel[pw] = 0;
				if (pe > 0) newLabel[pe] = 0;
			}
		}

		// north and south faces
		final int lastRow = w * (h - 1);
		for (int z = 0; z < d; z++) {
			for (int x = 0; x < w; x++) {
				final int pn = particleLabels[z][x];
				final int ps = particleLabels[z][lastRow + x];
				if (pn > 0) newLabel[pn] = 0;
				if (ps > 0) newLabel[ps] = 0;
			}
		}

		// replace labels
		final int wh = w * h;
		for (int z = 0; z < d; z++) {
			for (int i = 0; i < wh; i++) {
				final int p = particleLabels[z][i];
				final int nL = newLabel[p];
				if (nL == 0) {
					particleLabels[z][i] = 0;
					workArray[z][i] = 0;
				}
			}
		}

	}

	/**
	 * Remove particles outside user-specified volume thresholds
	 *
	 * @param imp ImagePlus, used for calibration
	 * @param workArray binary foreground and background information
	 * @param particleLabels Packed 3D array of particle labels
	 * @param minVol minimum (inclusive) particle volume
	 * @param maxVol maximum (inclusive) particle volume
	 * @param phase phase we are interested in
	 */
	static void filterParticles(final ImagePlus imp, final byte[][] workArray,
		final int[][] particleLabels, final double minVol, final double maxVol,
		final int phase)
	{
		if (minVol == 0 && maxVol == Double.POSITIVE_INFINITY) return;
		final int d = imp.getImageStackSize();
		final int wh = workArray[0].length;
		final long[] particleSizes = getParticleSizes(particleLabels);
		final double[] particleVolumes = getVolumes(imp, particleSizes);
		final byte flip;
		if (phase == FORE) {
			flip = 0;
		}
		else {
			flip = (byte) 255;
		}
		for (int z = 0; z < d; z++) {
			for (int i = 0; i < wh; i++) {
				final int p = particleLabels[z][i];
				final double v = particleVolumes[p];
				if (v < minVol || v > maxVol) {
					workArray[z][i] = flip;
					particleLabels[z][i] = 0;
				}
			}
		}
	}
}
