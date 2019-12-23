/*
BSD 2-Clause License
Copyright (c) 2018, Michael Doube, Richard Domander, Alessandro Felder
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

package org.bonej.plugins;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.TextField;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

import org.bonej.geometry.Ellipsoid;
import org.bonej.geometry.FitEllipsoid;
import org.bonej.menuWrappers.LocalThickness;
import org.bonej.util.DialogModifier;
import org.bonej.util.ImageCheck;
import org.bonej.util.Multithreader;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import customnode.CustomPointMesh;
import customnode.CustomTriangleMesh;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij3d.Image3DUniverse;
import marchingcubes.MCTriangulator;

/**
 * <p>
 * This class implements mutithreaded linear O(n) 3D particle
 * identification and shape analysis. It is a two-pass connected components labelling
 * algorithm, which uses reduction of a neighbour network to generate a lut.
 * Processing time increases linearly with number of pixels.
 * </p>
 *
 * @author Michael Doube
 */
public class ParticleCounter implements PlugIn, DialogListener {

//	/** Foreground value */
//	static final int FORE = -1;
//	/** Background value */
//	static final int BACK = 0;
	/** Surface colour style */
	private static final int GRADIENT = 0;
	private static final int SPLIT = 1;
	private static final int ORIENTATION = 2;
//	/** 2^23 - greatest integer that can be represented precisely by a float */
//	private static final int MAX_LABEL = 8388608; 
	/** String representation of current analysis phase for GUI display*/
	private String sPhase = "";

//TODO -- Run method & GUI (userland class)
	
	@Override
	public boolean dialogItemChanged(final GenericDialog gd, final AWTEvent e) {
		if (DialogModifier.hasInvalidNumber(gd.getNumericFields())) return false;
		final List<?> choices = gd.getChoices();
		final List<?> checkboxes = gd.getCheckboxes();
		final List<?> numbers = gd.getNumericFields();

		// link moments and ellipsoid choice to unit vector choice
		final Checkbox momBox = (Checkbox) checkboxes.get(4);
		final Checkbox elBox = (Checkbox) checkboxes.get(8);
		final Checkbox vvvBox = (Checkbox) checkboxes.get(9);
		vvvBox.setEnabled(elBox.getState() || momBox.getState());
		// link show stack 3d to volume resampling
		final Checkbox box = (Checkbox) checkboxes.get(17);
		final TextField numb = (TextField) numbers.get(4);
		numb.setEnabled(box.getState());
		// link show surfaces, gradient choice and split value
		final Checkbox surfbox = (Checkbox) checkboxes.get(13);
		final Choice col = (Choice) choices.get(0);
		final TextField split = (TextField) numbers.get(3);
		col.setEnabled(surfbox.getState());
		split.setEnabled(surfbox.getState() && col.getSelectedIndex() == 1);
		DialogModifier.registerMacroValues(gd, gd.getComponents());
		return true;
	}

	@Override
	public void run(final String arg) {
		final ImagePlus imp = IJ.getImage();
		if (null == imp) {
			IJ.noImage();
			return;
		}
		if (!ImageCheck.isBinary(imp)) {
			IJ.error("Binary image required");
			return;
		}
		final Calibration cal = imp.getCalibration();
		final String units = cal.getUnits();
		final GenericDialog gd = new GenericDialog("Setup");
		final String[] headers = { "Measurement Options", " " };
		final String[] labels = new String[10];
		final boolean[] defaultValues = new boolean[10];
		labels[0] = "Exclude on sides";
		defaultValues[0] = false;
		labels[1] = "Surface_area";
		defaultValues[1] = true;
		labels[2] = "Feret diameter";
		defaultValues[2] = false;
		labels[3] = "Enclosed_volume";
		defaultValues[3] = true;
		labels[4] = "Moments of inertia";
		defaultValues[4] = true;
		labels[5] = "Euler characteristic";
		defaultValues[5] = true;
		labels[6] = "Thickness";
		defaultValues[6] = true;
		labels[7] = "Mask thickness map";
		defaultValues[7] = false;
		labels[8] = "Ellipsoids";
		defaultValues[8] = true;
		labels[9] = "Record unit vectors";
		defaultValues[9] = false;
		gd.addCheckboxGroup(5, 2, labels, defaultValues, headers);
		gd.addNumericField("Min Volume", 0, 3, 7, units + "³");
		gd.addNumericField("Max Volume", Double.POSITIVE_INFINITY, 3, 7, units +
			"³");
		gd.addNumericField("Surface_resampling", 2, 0);
		final String[] headers2 = { "Graphical Results", " " };
		final String[] labels2 = new String[9];
		final boolean[] defaultValues2 = new boolean[9];
		labels2[0] = "Show_particle stack";
		defaultValues2[0] = true;
		labels2[1] = "Show_size stack";
		defaultValues2[1] = false;
		labels2[2] = "Show_thickness stack";
		defaultValues2[2] = false;
		labels2[3] = "Show_surfaces (3D)";
		defaultValues2[3] = true;
		labels2[4] = "Show_centroids (3D)";
		defaultValues2[4] = true;
		labels2[5] = "Show_axes (3D)";
		defaultValues2[5] = true;
		labels2[6] = "Show_ellipsoids (3D)";
		defaultValues2[6] = true;
		labels2[7] = "Show_stack (3D)";
		defaultValues2[7] = true;
		labels2[8] = "Draw_ellipsoids";
		defaultValues2[8] = false;
		gd.addCheckboxGroup(5, 2, labels2, defaultValues2, headers2);
		final String[] items = { "Gradient", "Split", "Orientation"};
		gd.addChoice("Surface colours", items, items[0]);
		gd.addNumericField("Split value", 0, 3, 7, units + "³");
		gd.addNumericField("Volume_resampling", 2, 0);

		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		final double minVol = gd.getNextNumber();
		final double maxVol = gd.getNextNumber();
		final boolean doExclude = gd.getNextBoolean();
		final boolean doSurfaceArea = gd.getNextBoolean();
		final boolean doFeret = gd.getNextBoolean();
		final boolean doSurfaceVolume = gd.getNextBoolean();
		final int resampling = (int) Math.floor(gd.getNextNumber());
		final boolean doMoments = gd.getNextBoolean();
		final boolean doEulerCharacters = gd.getNextBoolean();
		final boolean doThickness = gd.getNextBoolean();
		final boolean doMask = gd.getNextBoolean();
		final boolean doEllipsoids = gd.getNextBoolean();
		final boolean doVerboseUnitVectors = gd.getNextBoolean();
		final boolean doParticleImage = gd.getNextBoolean();
		final boolean doParticleSizeImage = gd.getNextBoolean();
		final boolean doThickImage = gd.getNextBoolean();
		final boolean doSurfaceImage = gd.getNextBoolean();
		final int colourMode = gd.getNextChoiceIndex();
		final double splitValue = gd.getNextNumber();
		final boolean doCentroidImage = gd.getNextBoolean();
		final boolean doAxesImage = gd.getNextBoolean();
		final boolean doEllipsoidImage = gd.getNextBoolean();
		final boolean do3DOriginal = gd.getNextBoolean();
		final boolean doEllipsoidStack = gd.getNextBoolean();
		final int origResampling = (int) Math.floor(gd.getNextNumber());

		// get the particles and do the analysis
		final long start = System.nanoTime();
		ConnectedComponents connector = new ConnectedComponents();
		final Object[] result = getParticles(connector, imp, minVol, maxVol,
				ConnectedComponents.FORE, doExclude);
		// calculate particle labelling time in ms
		final long time = (System.nanoTime() - start) / 1000000;
		IJ.log("Particle labelling finished in " + time + " ms");
		
		//start of analysis
		final int[][] particleLabels = (int[][]) result[1];
		final long[] particleSizes = (long[]) result[2];
		final int nParticles = connector.getNParticles();

		final double[] volumes = getVolumes(imp, particleSizes);
		final double[][] centroids = getCentroids(imp, particleLabels,
			particleSizes);
		final int[][] limits = getParticleLimits(imp, particleLabels, connector.getNParticles());

		// set up resources for analysis
		ArrayList<List<Point3f>> surfacePoints = new ArrayList<>();
		if (doSurfaceArea || doSurfaceVolume || doSurfaceImage || doEllipsoids ||
			doFeret || doEllipsoidStack)
		{
			surfacePoints = getSurfacePoints(imp, particleLabels, limits, resampling, connector.getNParticles());
		}
		EigenvalueDecomposition[] eigens = new EigenvalueDecomposition[nParticles];
		if (doMoments || doAxesImage || colourMode == ORIENTATION) {
			eigens = getEigens(imp, particleLabels, centroids, connector.getNParticles());
		}
		// calculate dimensions
		double[] surfaceAreas = new double[nParticles];
		if (doSurfaceArea) {
			surfaceAreas = getSurfaceAreas(surfacePoints);
		}
		double[] ferets = new double[nParticles];
		if (doFeret) {
			ferets = getFerets(surfacePoints);
		}
		double[] surfaceVolumes = new double[nParticles];
		if (doSurfaceVolume) {
			surfaceVolumes = getSurfaceVolume(surfacePoints);
		}
		double[][] eulerCharacters = new double[nParticles][3];
		if (doEulerCharacters) {
			eulerCharacters = getEulerCharacter(imp, particleLabels, limits, connector.getNParticles());
		}
		double[][] thick = new double[nParticles][2];
		if (doThickness) {
			final LocalThickness th = new LocalThickness();
			final ImagePlus thickImp = th.getLocalThickness(imp, false, doMask);
			thick = getMeanStdDev(thickImp, particleLabels, particleSizes);
			if (doThickImage) {
				double max = 0;
				for (int i = 1; i < nParticles; i++) {
					max = Math.max(max, thick[i][2]);
				}
				thickImp.getProcessor().setMinAndMax(0, max);
				thickImp.setTitle(imp.getShortTitle() + "_thickness");
				thickImp.show();
				thickImp.setSlice(1);
				IJ.run("Fire");
			}
		}
		Object[][] ellipsoids = new Object[nParticles][10];
		if (doEllipsoids || doEllipsoidImage || doEllipsoidStack) {
			ellipsoids = getEllipsoids(surfacePoints);
		}

		// Show numerical results
		final ResultsTable rt = new ResultsTable();
		for (int i = 1; i < volumes.length; i++) {
			if (volumes[i] > 0) {
				rt.incrementCounter();
				rt.addLabel(imp.getTitle());
				rt.addValue("ID", i);
				rt.addValue("Vol. (" + units + "³)", volumes[i]);
				rt.addValue("x Cent (" + units + ")", centroids[i][0]);
				rt.addValue("y Cent (" + units + ")", centroids[i][1]);
				rt.addValue("z Cent (" + units + ")", centroids[i][2]);
				if (doSurfaceArea) {
					rt.addValue("SA (" + units + "²)", surfaceAreas[i]);
				}
				if (doFeret) {
					rt.addValue("Feret (" + units + ")", ferets[i]);
				}
				if (doSurfaceVolume) {
					rt.addValue("Encl. Vol. (" + units + "³)", surfaceVolumes[i]);
				}
				if (doMoments) {
					final EigenvalueDecomposition E = eigens[i];
					rt.addValue("I1", E.getD().get(2, 2));
					rt.addValue("I2", E.getD().get(1, 1));
					rt.addValue("I3", E.getD().get(0, 0));
					rt.addValue("vX", E.getV().get(0, 0));
					rt.addValue("vY", E.getV().get(1, 0));
					rt.addValue("vZ", E.getV().get(2, 0));
					if (doVerboseUnitVectors) {
						rt.addValue("vX1", E.getV().get(0, 1));
						rt.addValue("vY1", E.getV().get(1, 1));
						rt.addValue("vZ1", E.getV().get(2, 1));
						rt.addValue("vX2", E.getV().get(0, 2));
						rt.addValue("vY2", E.getV().get(1, 2));
						rt.addValue("vZ2", E.getV().get(2, 2));
					}
				}
				if (doEulerCharacters) {
					rt.addValue("Euler (χ)", eulerCharacters[i][0]);
					rt.addValue("Holes (β1)", eulerCharacters[i][1]);
					rt.addValue("Cavities (β2)", eulerCharacters[i][2]);
				}
				if (doThickness) {
					rt.addValue("Thickness (" + units + ")", thick[i][0]);
					rt.addValue("SD Thickness (" + units + ")", thick[i][1]);
					rt.addValue("Max Thickness (" + units + ")", thick[i][2]);
				}
				if (doEllipsoids) {
					final double[] rad;
					final double[][] unitV;
					if (ellipsoids[i] == null) {
						rad = new double[] { Double.NaN, Double.NaN, Double.NaN };
						unitV = new double[][] { { Double.NaN, Double.NaN, Double.NaN }, {
							Double.NaN, Double.NaN, Double.NaN }, { Double.NaN, Double.NaN,
								Double.NaN } };
					}
					else {
						final Object[] el = ellipsoids[i];
						rad = (double[]) el[1];
						unitV = (double[][]) el[2];
					}
					rt.addValue("Major radius (" + units + ")", rad[0]);
					rt.addValue("Int. radius (" + units + ")", rad[1]);
					rt.addValue("Minor radius (" + units + ")", rad[2]);
					if (doVerboseUnitVectors) {
						rt.addValue("V00", unitV[0][0]);
						rt.addValue("V01", unitV[0][1]);
						rt.addValue("V02", unitV[0][2]);
						rt.addValue("V10", unitV[1][0]);
						rt.addValue("V11", unitV[1][1]);
						rt.addValue("V12", unitV[1][2]);
						rt.addValue("V20", unitV[2][0]);
						rt.addValue("V21", unitV[2][1]);
						rt.addValue("V22", unitV[2][2]);
					}
				}
				rt.updateResults();
			}
		}
		rt.show("Results");

		// Show resulting image stacks
		if (doParticleImage) {
			displayParticleLabels(particleLabels, imp).show();
			IJ.run("3-3-2 RGB");
		}
		if (doParticleSizeImage) {
			displayParticleValues(imp, particleLabels, volumes).show();
			IJ.run("Fire");
		}
		if (doEllipsoidStack) {
			displayParticleEllipsoids(imp, ellipsoids, "Ellipsoids").show();
		}

		// show 3D renderings
		if (doSurfaceImage || doCentroidImage || doAxesImage || do3DOriginal ||
			doEllipsoidImage)
		{

			final Image3DUniverse univ = new Image3DUniverse();
			if (doSurfaceImage) {
				displayParticleSurfaces(univ, surfacePoints, colourMode, volumes,
					splitValue, eigens);
			}
			if (doCentroidImage) {
				displayCentroids(centroids, univ);
			}
			if (doAxesImage) {
				displayPrincipalAxes(univ, eigens, centroids, particleSizes);
			}
			if (doEllipsoidImage) {
				displayEllipsoids(ellipsoids, univ);
			}
			if (do3DOriginal) {
				display3DOriginal(imp, origResampling, univ);
			}
			univ.show();
		}
		IJ.showProgress(1.0);
		IJ.showStatus("Particle Analysis Complete");
		UsageReporter.reportEvent(this).send();
	}
	
	//TODO--------Connected components labelling
	
	/**
	 * Get particles, particle labels and particle sizes from a 3D ImagePlus
	 *
	 * @param imp Binary input image
	 * @param minVol minimum volume particle to include
	 * @param maxVol maximum volume particle to include
	 * @param phase foreground or background (FORE or BACK)
	 * @param doExclude if true, remove particles touching sides of the stack
	 * @return Object[] {byte[][], int[][]} containing a binary workArray and
	 *         particle labels.
	 */
	private Object[] getParticles(ConnectedComponents connector, final ImagePlus imp,
		final double minVol, final double maxVol, final int phase,
		final boolean doExclude)
	{
		final byte[][] workArray = makeWorkArray(imp);
		return getParticles(connector, imp, workArray, minVol, maxVol, phase,
			doExclude);
	}

	/**
	 * Get particles, particle labels and particle sizes from a 3D ImagePlus
	 *  
	 * @param imp
	 * @param phase
	 * @return
	 */
	Object[] getParticles(ConnectedComponents connector, final ImagePlus imp, final int phase)
	{
		final byte[][] workArray = makeWorkArray(imp);
		return getParticles(connector, imp, workArray, 0.0,
			Double.POSITIVE_INFINITY, phase, false);
	}

	/**
	 * Get particles, particle labels and particle sizes from a 3D ImagePlus
	 * 
	 * @param imp
	 * @param workArray
	 * @param phase
	 * @return
	 */
	Object[] getParticles(ConnectedComponents connector, 
			final ImagePlus imp, final byte[][] workArray, final int phase)
	{
		return getParticles(connector, imp, workArray, 0.0,
			Double.POSITIVE_INFINITY, phase, false);
	}
	
	/**
	 * Get particles, particle labels and sizes from a workArray using an
	 * ImagePlus for scale information
	 *
	 * @param imp input binary image
	 * @param workArray work array
	 * @param minVol minimum volume particle to include
	 * @param maxVol maximum volume particle to include
	 * @param phase FORE or BACK for foreground or background respectively
	 * @param doExclude exclude particles touching the edges.
	 * @return Object[] array containing a binary workArray, particle labels and
	 *         particle sizes
	 */
	private Object[] getParticles(ConnectedComponents connector, 
			final ImagePlus imp, final byte[][] workArray,
		final double minVol, final double maxVol,
		final int phase, final boolean doExclude)
	{
		if (phase == ConnectedComponents.FORE) {
			sPhase = "foreground";
		}
		else if (phase == ConnectedComponents.BACK) {
			sPhase = "background";
		}
		else {
			throw new IllegalArgumentException();
		}
		
		//first pass through whole stack
		final int[][] particleLabels = connector.firstIDAttribution(imp, workArray, phase);

		final int nParticles = connector.getNParticles();
		
		filterParticles(imp, workArray, particleLabels, minVol, maxVol, phase, nParticles);
		
		if (doExclude) excludeOnEdges(imp, particleLabels, workArray, nParticles);

		final long[] particleSizes = getParticleSizes(particleLabels, connector.getNParticles());
		return new Object[] { workArray, particleLabels, particleSizes };
	}
	
	/**
	 * Create a work array
	 *
	 * @param imp an image.
	 * @return byte[][] work array
	 */
	private static byte[][] makeWorkArray(final ImagePlus imp) {
		final int s = imp.getStackSize();
		final int p = imp.getWidth() * imp.getHeight();
		final byte[][] workArray = new byte[s][p];
		final ImageStack stack = imp.getStack();

		AtomicInteger ai = new AtomicInteger(0);

		final Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(() -> {
				for (int z = ai.getAndIncrement(); z < s; z = ai.getAndIncrement()) {
					final ImageProcessor ip = stack.getProcessor(z + 1);
					for (int i = 0; i < p; i++) {
						workArray[z][i] = (byte) ip.get(i);
					}
				}
			});
		}
		Multithreader.startAndJoin(threads);
		return workArray;
	}
	
	
	
//TODO------------ANALYSIS---------REMOVE-TO-OWN-CLASS-------------------------
	
	/**
	 * Get the sizes of all the particles as a voxel count
	 *
	 * @param particleLabels particles in the image.
	 * @return particleSizes sizes of the particles.
	 */
	long[] getParticleSizes(final int[][] particleLabels, final int nParticles) {
		IJ.showStatus("Getting " + sPhase + " particle sizes");
		final int d = particleLabels.length;
		final int wh = particleLabels[0].length;
		
		//make a list of all the particle sizes with 
		//index = particle value
		//need to handle the ID offset for the chunks
		AtomicInteger an = new AtomicInteger(0);
		final long[][] partSizes = new long[d][];
		
		final Thread[] threads = Multithreader.newThreads();
		for (int thread = 0; thread < threads.length; thread++) {
			threads[thread] = new Thread(() -> {
				for (int z = an.getAndIncrement(); z < d; z = an.getAndIncrement()) {
					final long[] particleSizes = new long[nParticles];
					final int[] slice = particleLabels[z];
					for (int i = 0; i < wh; i++) {
						final int label = slice[i];
						//hack to avoid AIOOB
//						if (label <= maxParticle)
							particleSizes[label]++;
					}
					partSizes[z] = particleSizes;
				}
			});
		}
		Multithreader.startAndJoin(threads);
		
		final long[] particleSizes = new long[nParticles];
		for (int i = 0; i < nParticles; i++) {
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

	private static double[] getSurfaceAreas(
		final Collection<List<Point3f>> surfacePoints)
	{
		return surfacePoints.stream().mapToDouble(ParticleCounter::getSurfaceArea)
			.toArray();
	}

	private static ArrayList<List<Point3f>> getSurfacePoints(final ImagePlus imp,
		final int[][] particleLabels, final int[][] limits, final int resampling, final int nParticles)
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

	private static double[] getSurfaceVolume(
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

	private static double[] getVolumes(final ImagePlus imp,
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

	private static boolean isRadiiValid(final double[] radii) {
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

	private static void display3DOriginal(final ImagePlus imp,
		final int resampling, final Image3DUniverse univ)
	{
		final Color3f colour = new Color3f(1.0f, 1.0f, 1.0f);
		final boolean[] channels = { true, true, true };
		try {
			univ.addVoltex(imp, colour, imp.getTitle(), 0, channels, resampling)
				.setLocked(true);
		}
		catch (final NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
		}
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
		IJ.log("Copying for analysis particle ID="+p+" with xMin="+xMin+", xMax="+xMax
				+ " yMin="+yMin+", yMax="+yMax+", zMin="+zMin+", zMax="+zMax);
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
	private static double[][] getCentroids(final ImagePlus imp,
		final int[][] particleLabels, final long[] particleSizes)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int nParticles = particleSizes.length;
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

	private static EigenvalueDecomposition[] getEigens(final ImagePlus imp,
		final int[][] particleLabels, final double[][] centroids, final int nParticles)
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

	private static Object[][] getEllipsoids(
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
	private double[][] getEulerCharacter(final ImagePlus imp,
		final int[][] particleLabels, final int[][] limits, final int nParticles)
	{
		final Connectivity con = new Connectivity();
		final ConnectedComponents conComp = new ConnectedComponents();
		final double[][] eulerCharacters = new double[nParticles][3];
		for (int p = 1; p < nParticles; p++) {
			final ImagePlus particleImp = getBinaryParticle(p, imp, particleLabels,
				limits, 1);
			final double euler = con.getSumEuler(particleImp);
			final double cavities = getNCavities(conComp, particleImp);
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
	private static double[] getFerets(
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
	private static double[][] getMeanStdDev(final ImagePlus imp,
		final int[][] particleLabels, final long[] particleSizes)
	{
		final int d = imp.getImageStackSize();
		final int wh = imp.getWidth() * imp.getHeight();
		final ImageStack stack = imp.getImageStack();
		final int nParticles = particleSizes.length;
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

	private int getNCavities(ConnectedComponents connector, final ImagePlus imp) {
		final Object[] result = getParticles(connector, imp, ConnectedComponents.BACK);
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
	private static int[][] getParticleLimits(final ImagePlus imp,
		final int[][] particleLabels, final int nParticles)
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
	private void excludeOnEdges(final ImagePlus imp, final int[][] particleLabels,
		final byte[][] workArray, final int nParticles)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final long[] particleSizes = getParticleSizes(particleLabels, nParticles);
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
	private void filterParticles(final ImagePlus imp, final byte[][] workArray,
		final int[][] particleLabels, final double minVol, final double maxVol,
		final int phase, final int nParticles)
	{
		if (minVol == 0 && maxVol == Double.POSITIVE_INFINITY) return;
		final int d = imp.getImageStackSize();
		final int wh = workArray[0].length;
		final long[] particleSizes = getParticleSizes(particleLabels, nParticles);
		final double[] particleVolumes = getVolumes(imp, particleSizes);
		final byte flip;
		if (phase == ConnectedComponents.FORE) {
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
//TODO------------DISPLAY---------REMOVE-TO-OWN-CLASS-------------------------
	
	/**
	 * Generate a colour based on the inertia tensor's eigenvector
	 * 
	 * Colour is from the HSB colour wheel scaled by 0.5 to fit
	 * into pi radians (rather than the 2 pi it normally occupies),
	 * so that red is at 0, pi and 2pi radians.
	 * 
	 * Colour is mapped to the axis-angle representation of the tensor
	 * so hue varies as a function of second axis rotation around the
	 * first.
	 * 
	 * @param eigen Eigenvalue decomposition of the particle
	 * @return Colour scaling in red for axis and green for angle
	 */
	private static Color3f colourFromEigenVector(EigenvalueDecomposition eigen)
	{
		final Matrix rotation = eigen.getV();
		
		//deflection of long axis from image z axis, 0 - pi radians
		final double angle = Math.acos(-Math.abs(rotation.get(2, 0)));
		
	  final float hue = (float)(angle / Math.PI);
		final float saturation = 1.0f;
		final float brightness = 1.0f;
		
		final int rgb = Color.HSBtoRGB(hue, saturation, brightness);
		final Color color = new Color(rgb);
		float red = (float)(color.getRed()/255d);
    float green = (float)(color.getGreen()/255d);
    float blue = (float)(color.getBlue()/255d);
	
		return new Color3f(red, green, blue);
	}

	/**
	 * Create an image showing some particle measurement
	 *
	 * @param imp an image.
	 * @param particleLabels the particles in the image.
	 * @param values list of values whose array indices correspond to
	 *          particlelabels
	 * @return ImagePlus with particle labels substituted with some value
	 */
	private static ImagePlus displayParticleValues(final ImagePlus imp,
		final int[][] particleLabels, final double[] values)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int wh = w * h;
		final float[][] pL = new float[d][wh];
		values[0] = 0; // don't colour the background
		final ImageStack stack = new ImageStack(w, h);
		for (int z = 0; z < d; z++) {
			for (int i = 0; i < wh; i++) {
				final int p = particleLabels[z][i];
				pL[z][i] = (float) values[p];
			}
			stack.addSlice(imp.getImageStack().getSliceLabel(z + 1), pL[z]);
		}
		final double max = Arrays.stream(values).max().orElse(0.0);
		final ImagePlus impOut = new ImagePlus(imp.getShortTitle() + "_" + "volume",
			stack);
		impOut.setCalibration(imp.getCalibration());
		impOut.getProcessor().setMinAndMax(0, max);
		return impOut;
	}

	/**
	 * 
	 * @param imp
	 * @param ellipsoids
	 * @param title
	 * @return ImagePlus containing particles drawn as best-fit solid ellipsoids
	 */
	private ImagePlus displayParticleEllipsoids(final ImagePlus imp, final Object[][] ellipsoids,
		final String title) {
	final int w = imp.getWidth();
	final int h = imp.getHeight();
	final int d = imp.getImageStackSize();
	
	Calibration cal = imp.getCalibration();
	final double pW = cal.pixelWidth;
	final double pH = cal.pixelHeight;
	final double pD = cal.pixelDepth;
	
	//set up a work array
	final ByteProcessor[] bps = new ByteProcessor[d];
	for (int z = 0; z < d; z++) {
		bps[z] = new ByteProcessor(w, h);
	}
	
	final int n = ellipsoids.length;
	for (int i = 0; i < n; i++) {
		IJ.showStatus("Drawing ellipsoid stack...");
		IJ.showProgress(i, n);
		Ellipsoid ellipsoid;
		try	{
			ellipsoid = new Ellipsoid(ellipsoids[i]);
		} catch (Exception e) {
			continue;
		}
		
		//ellipsoid is in calibrated real-world units
		final double[] box = ellipsoid.getAxisAlignedBoundingBox();
		
		//decalibrate to pixels
		final int xMin = clamp((int) Math.floor(box[0] / pW), 0, w-1);
		final int xMax = clamp((int) Math.floor(box[1] / pW), 0, w-1);
		final int yMin = clamp((int) Math.floor(box[2] / pH), 0, h-1);
		final int yMax = clamp((int) Math.floor(box[3] / pH), 0, h-1);
		final int zMin = clamp((int) Math.floor(box[4] / pD), 0, d-1);
		final int zMax = clamp((int) Math.floor(box[5] / pD), 0, d-1);
		
		//set the ellipsoid-contained pixels to foreground
		for (int z = zMin; z <= zMax; z++) {
			for (int y = yMin; y <= yMax; y++) {
				for (int x = xMin; x <= xMax; x++ ) {
					if (ellipsoid.contains(x * pW, y * pH, z * pD)){
						bps[z].set(x, y, 255);
					}
				}
			}
		}		
	}
	
	ImageStack stack = new ImageStack(w, h);
	for (ByteProcessor bp : bps)
		stack.addSlice(bp);
	
	final ImagePlus impOut = new ImagePlus(imp.getShortTitle() + "_" + title, stack);
	impOut.setCalibration(cal);
	return impOut;
  }
	
	private int clamp(int value, int min, int max) {
		if (value < min)
			return min;
		if (value > max)
			return max;
		return value;
	}


	/**
	 * Draws 3 orthogonal axes defined by the centroid, unitvector and axis
	 * length.
	 *
	 * @param univ the universe where axes are drawn.
	 * @param centroid centroid of a particle.
	 * @param unitVector orientation of the particle.
	 * @param lengths lengths of the axes.
	 * @param green green component of the axes' color.
	 * @param title text shown by the axes.
	 */
	private static void displayAxes(final Image3DUniverse univ,
		final double[] centroid, final double[][] unitVector,
		final double[] lengths, final float green, final String title)
	{
		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];
		final double eVec1x = unitVector[0][0];
		final double eVec1y = unitVector[1][0];
		final double eVec1z = unitVector[2][0];
		final double eVec2x = unitVector[0][1];
		final double eVec2y = unitVector[1][1];
		final double eVec2z = unitVector[2][1];
		final double eVec3x = unitVector[0][2];
		final double eVec3y = unitVector[1][2];
		final double eVec3z = unitVector[2][2];
		final double l1 = lengths[0];
		final double l2 = lengths[1];
		final double l3 = lengths[2];

		final List<Point3f> mesh = new ArrayList<>();
		final Point3f start1 = new Point3f();
		start1.x = (float) (cX - eVec1x * l1);
		start1.y = (float) (cY - eVec1y * l1);
		start1.z = (float) (cZ - eVec1z * l1);
		mesh.add(start1);

		final Point3f end1 = new Point3f();
		end1.x = (float) (cX + eVec1x * l1);
		end1.y = (float) (cY + eVec1y * l1);
		end1.z = (float) (cZ + eVec1z * l1);
		mesh.add(end1);

		final Point3f start2 = new Point3f();
		start2.x = (float) (cX - eVec2x * l2);
		start2.y = (float) (cY - eVec2y * l2);
		start2.z = (float) (cZ - eVec2z * l2);
		mesh.add(start2);

		final Point3f end2 = new Point3f();
		end2.x = (float) (cX + eVec2x * l2);
		end2.y = (float) (cY + eVec2y * l2);
		end2.z = (float) (cZ + eVec2z * l2);
		mesh.add(end2);

		final Point3f start3 = new Point3f();
		start3.x = (float) (cX - eVec3x * l3);
		start3.y = (float) (cY - eVec3y * l3);
		start3.z = (float) (cZ - eVec3z * l3);
		mesh.add(start3);

		final Point3f end3 = new Point3f();
		end3.x = (float) (cX + eVec3x * l3);
		end3.y = (float) (cY + eVec3y * l3);
		end3.z = (float) (cZ + eVec3z * l3);
		mesh.add(end3);

		final Color3f aColour = new Color3f(1.0f, green, 0.0f);
		try {
			univ.addLineMesh(mesh, aColour, title, false).setLocked(true);
		}
		catch (final NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
		}
	}

	/**
	 * Draw the particle centroids in a 3D viewer
	 *
	 * @param centroids [n][3] centroids of particles.
	 * @param univ universe where the centroids are displayed.
	 */
	private static void displayCentroids(final double[][] centroids,
		final Image3DUniverse univ)
	{
		final int nCentroids = centroids.length;
		for (int p = 1; p < nCentroids; p++) {
			IJ.showStatus("Rendering centroids...");
			IJ.showProgress(p, nCentroids);
			final Point3f centroid = new Point3f();
			centroid.x = (float) centroids[p][0];
			centroid.y = (float) centroids[p][1];
			centroid.z = (float) centroids[p][2];
			final List<Point3f> point = new ArrayList<>();
			point.add(centroid);
			final CustomPointMesh mesh = new CustomPointMesh(point);
			mesh.setPointSize(5.0f);
			final float red = 0.0f;
			final float green = 0.5f * p / nCentroids;
			final float blue = 1.0f;
			final Color3f cColour = new Color3f(red, green, blue);
			mesh.setColor(cColour);
			try {
				univ.addCustomMesh(mesh, "Centroid " + p).setLocked(true);
			}
			catch (final NullPointerException npe) {
				IJ.log("3D Viewer was closed before rendering completed.");
				return;
			}
		}
	}

	private static void displayEllipsoids(final Object[][] ellipsoids,
		final Image3DUniverse univ)
	{
		final int nEllipsoids = ellipsoids.length;
		for (int el = 1; el < nEllipsoids; el++) {
			IJ.showStatus("Rendering ellipsoids...");
			IJ.showProgress(el, nEllipsoids);
			if (ellipsoids[el] == null) {
				continue;
			}
			final double[] radii = (double[]) ellipsoids[el][1];
			if (!isRadiiValid(radii)) {
				continue;
			}
			final double[] centre = (double[]) ellipsoids[el][0];
			final double[][] eV = (double[][]) ellipsoids[el][2];
			final double a = radii[0]; // longest
			final double b = radii[1]; // middle
			final double c = radii[2]; // shortest
			if (a < b || b < c || a < c) {
				IJ.log("Error: Bad ellipsoid radius ordering! Surface: " + el);
			}
			final double[][] ellipsoid = FitEllipsoid.testEllipsoid(a, b, c, 0, 0, 0,
				0, 0, 1000, false);
			final int nPoints = ellipsoid.length;
			// rotate points by eigenvector matrix
			// and add transformation for centre
			for (int p = 0; p < nPoints; p++) {
				final double x = ellipsoid[p][0];
				final double y = ellipsoid[p][1];
				final double z = ellipsoid[p][2];
				ellipsoid[p][0] = x * eV[0][0] + y * eV[0][1] + z * eV[0][2] +
					centre[0];
				ellipsoid[p][1] = x * eV[1][0] + y * eV[1][1] + z * eV[1][2] +
					centre[1];
				ellipsoid[p][2] = x * eV[2][0] + y * eV[2][1] + z * eV[2][2] +
					centre[2];
			}

			final List<Point3f> points = new ArrayList<>();
			for (final double[] anEllipsoid : ellipsoid) {
				final Point3f e = new Point3f();
				e.x = (float) anEllipsoid[0];
				e.y = (float) anEllipsoid[1];
				e.z = (float) anEllipsoid[2];
				points.add(e);
			}
			final CustomPointMesh mesh = new CustomPointMesh(points);
			mesh.setPointSize(1.0f);
			final float red = 0.0f;
			final float green = 0.5f;
			final float blue = 1.0f;
			final Color3f cColour = new Color3f(red, green, blue);
			mesh.setColor(cColour);
			try {
				univ.addCustomMesh(mesh, "Ellipsoid " + el).setLocked(true);
			}
			catch (final NullPointerException npe) {
				IJ.log("3D Viewer was closed before rendering completed.");
				return;
			}
			// Add some axes
			displayAxes(univ, centre, eV, radii, 1.0f, "Ellipsoid Axes " + el);
		}
	}

	/**
	 * Display the particle labels as an ImagePlus
	 *
	 * @param particleLabels particles labelled in the original image.
	 * @param imp original image, used for image dimensions, calibration and
	 *          titles
	 * @return an image of the particles.
	 */
	private static ImagePlus displayParticleLabels(final int[][] particleLabels,
		final ImagePlus imp)
	{
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();
		final int wh = w * h;
		final ImageStack stack = new ImageStack(w, h);
		double max = 0;
		for (int z = 0; z < d; z++) {
			final float[] slicePixels = new float[wh];
			for (int i = 0; i < wh; i++) {
				slicePixels[i] = particleLabels[z][i];
				max = Math.max(max, slicePixels[i]);
			}
			stack.addSlice(imp.getImageStack().getSliceLabel(z + 1), slicePixels);
		}
		final ImagePlus impParticles = new ImagePlus(imp.getShortTitle() + "_parts",
			stack);
		impParticles.setCalibration(imp.getCalibration());
		impParticles.getProcessor().setMinAndMax(0, max);
		if (max > ConnectedComponents.MAX_LABEL) IJ.error("Warning",
			"More than 2^23 particles. " +
				"Particle label display is imprecise above this number due to int to float conversion.");
		return impParticles;
	}

	/**
	 * Draw the particle surfaces in a 3D viewer
	 *
	 * @param univ universe where the centroids are displayed.
	 * @param surfacePoints points of each particle.
	 */
	private static void displayParticleSurfaces(final Image3DUniverse univ,
		final Collection<List<Point3f>> surfacePoints, final int colourMode,
		final double[] volumes, final double splitValue,
		final EigenvalueDecomposition[] eigens)
	{
		int p = 0;
		final int nSurfaces = surfacePoints.size();
		for (final List<Point3f> surfacePoint : surfacePoints) {
			IJ.showStatus("Rendering surfaces...");
			IJ.showProgress(p, nSurfaces);
			if (p > 0 && !surfacePoint.isEmpty()) {
				Color3f pColour = new Color3f(0, 0, 0);
				if (colourMode == GRADIENT) {
					final float red = 1.0f - p / (float) nSurfaces;
					final float green = 1.0f - red;
					final float blue = p / (2.0f * nSurfaces);
					pColour = new Color3f(red, green, blue);
				}
				else if (colourMode == SPLIT) {
					if (volumes[p] > splitValue) {
						// red if over
						pColour = new Color3f(1.0f, 0.0f, 0.0f);
					}
					else {
						// yellow if under
						pColour = new Color3f(1.0f, 1.0f, 0.0f);
					}
				}
				else if (colourMode == ORIENTATION) {
					pColour = colourFromEigenVector(eigens[p]);
				}
				// Add the mesh
				try {
					univ.addTriangleMesh(surfacePoint, pColour, "Surface " + p).setLocked(
						true);
				}
				catch (final NullPointerException npe) {
					IJ.log("3D Viewer was closed before rendering completed.");
					return;
				}
			}
			p++;
		}
	}
	
	private static void displayPrincipalAxes(final Image3DUniverse univ,
		final EigenvalueDecomposition[] eigens, final double[][] centroids,
		long[] particleSizes)
	{
		final int nEigens = eigens.length;
				
		for (int p = 1; p < nEigens; p++) {
			IJ.showStatus("Rendering principal axes...");
			IJ.showProgress(p, nEigens);
			
			final long size = particleSizes[p];
			final Matrix eVec = eigens[p].getV();
			final Matrix eVal = eigens[p].getD();
			double[] lengths = new double[3];
			for (int i = 0; i < 3; i++) {
				lengths[i] = 2 * Math.sqrt(eVal.get(2 - i, 2 - i) / size);
			}
			displayAxes(univ, centroids[p], eVec.getArray(), lengths, 0.0f,
				"Principal Axes " + p);
		}
	}

}


