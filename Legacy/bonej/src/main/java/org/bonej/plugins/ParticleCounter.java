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
import java.awt.TextField;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import org.bonej.menuWrappers.LocalThickness;
import org.bonej.util.DialogModifier;
import org.bonej.util.ImageCheck;
import org.bonej.util.Multithreader;
import org.scijava.vecmath.Point3f;

import Jama.EigenvalueDecomposition;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij3d.Image3DUniverse;

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

	/** Foreground value */
	static final int FORE = -1;
	/** Background value */
	static final int BACK = 0;
	/** Surface colour style */
	private static final int GRADIENT = 0;
	private static final int SPLIT = 1;
	private static final int ORIENTATION = 2;
	/** 2^23 - greatest integer that can be represented precisely by a float */
	private static final int MAX_LABEL = 8388608; 
	private String sPhase = "";
	
	/** number of particle labels */
	private static int nParticles;

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
		final Object[] result = getParticles(imp, minVol, maxVol,	FORE, doExclude);
		// calculate particle labelling time in ms
		final long time = (System.nanoTime() - start) / 1000000;
		IJ.log("Particle labelling finished in " + time + " ms");
		
		//start of analysis
		final int[][] particleLabels = (int[][]) result[1];
		final long[] particleSizes = (long[]) result[2];
		nParticles = particleSizes.length;
		final double[] volumes = ParticleAnalysis.getVolumes(imp, particleSizes);
		final double[][] centroids = ParticleAnalysis.getCentroids(imp, particleLabels,
			particleSizes);
		final int[][] limits = ParticleAnalysis.getParticleLimits(imp, particleLabels);

		// set up resources for analysis
		ArrayList<List<Point3f>> surfacePoints = new ArrayList<>();
		if (doSurfaceArea || doSurfaceVolume || doSurfaceImage || doEllipsoids ||
			doFeret || doEllipsoidStack)
		{
			surfacePoints = ParticleAnalysis.getSurfacePoints(imp, particleLabels, limits, resampling);
		}
		EigenvalueDecomposition[] eigens = new EigenvalueDecomposition[nParticles];
		if (doMoments || doAxesImage || colourMode == ORIENTATION) {
			eigens = ParticleAnalysis.getEigens(imp, particleLabels, centroids);
		}
		// calculate dimensions
		double[] surfaceAreas = new double[nParticles];
		if (doSurfaceArea) {
			surfaceAreas = ParticleAnalysis.getSurfaceAreas(surfacePoints);
		}
		double[] ferets = new double[nParticles];
		if (doFeret) {
			ferets = ParticleAnalysis.getFerets(surfacePoints);
		}
		double[] surfaceVolumes = new double[nParticles];
		if (doSurfaceVolume) {
			surfaceVolumes = ParticleAnalysis.getSurfaceVolume(surfacePoints);
		}
		double[][] eulerCharacters = new double[nParticles][3];
		if (doEulerCharacters) {
			eulerCharacters = ParticleAnalysis.getEulerCharacter(imp, particleLabels, limits);
		}
		double[][] thick = new double[nParticles][2];
		if (doThickness) {
			final LocalThickness th = new LocalThickness();
			final ImagePlus thickImp = th.getLocalThickness(imp, false, doMask);
			thick = ParticleAnalysis.getMeanStdDev(thickImp, particleLabels, particleSizes);
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
			ellipsoids = ParticleAnalysis.getEllipsoids(surfacePoints);
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
			ParticleDisplay.displayParticleLabels(particleLabels, imp).show();
			IJ.run("3-3-2 RGB");
		}
		if (doParticleSizeImage) {
			ParticleDisplay.displayParticleValues(imp, particleLabels, volumes).show();
			IJ.run("Fire");
		}
		if (doEllipsoidStack) {
			ParticleDisplay.displayParticleEllipsoids(imp, ellipsoids, "Ellipsoids").show();
		}

		// show 3D renderings
		if (doSurfaceImage || doCentroidImage || doAxesImage || do3DOriginal ||
			doEllipsoidImage)
		{

			final Image3DUniverse univ = new Image3DUniverse();
			if (doSurfaceImage) {
				ParticleDisplay.displayParticleSurfaces(univ, surfacePoints, colourMode, volumes,
					splitValue, eigens);
			}
			if (doCentroidImage) {
				ParticleDisplay.displayCentroids(centroids, univ);
			}
			if (doAxesImage) {
				ParticleDisplay.displayPrincipalAxes(univ, eigens, centroids, particleSizes);
			}
			if (doEllipsoidImage) {
				ParticleDisplay.displayEllipsoids(ellipsoids, univ);
			}
			if (do3DOriginal) {
				ParticleDisplay.display3DOriginal(imp, origResampling, univ);
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
	private Object[] getParticles(final ImagePlus imp,
		final double minVol, final double maxVol, final int phase,
		final boolean doExclude)
	{
		final byte[][] workArray = makeWorkArray(imp);
		return getParticles(imp, workArray, minVol, maxVol, phase,
			doExclude);
	}

	/**
	 * Get particles, particle labels and particle sizes from a 3D ImagePlus
	 *  
	 * @param imp
	 * @param phase
	 * @return
	 */
	Object[] getParticles(final ImagePlus imp, final int phase)
	{
		final byte[][] workArray = makeWorkArray(imp);
		return getParticles(imp, workArray, 0.0,
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
	Object[] getParticles(final ImagePlus imp, final byte[][] workArray, final int phase)
	{
		return getParticles(imp, workArray, 0.0,
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
	private Object[] getParticles(final ImagePlus imp, final byte[][] workArray,
		final double minVol, final double maxVol,
		final int phase, final boolean doExclude)
	{
		if (phase == FORE) {
			sPhase = "foreground";
		}
		else if (phase == BACK) {
			sPhase = "background";
		}
		else {
			throw new IllegalArgumentException();
		}
		
		//first pass through whole stack
		final int[][] particleLabels = firstIDAttribution(imp, workArray, phase);

		ParticleAnalysis.filterParticles(imp, workArray, particleLabels, minVol, maxVol, phase);
		
		if (doExclude) ParticleAnalysis.excludeOnEdges(imp, particleLabels, workArray);

		final long[] particleSizes = ParticleAnalysis.getParticleSizes(particleLabels);
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
	
	/**
	 * Go through all pixels and assign initial particle label.
	 *
	 * @param imp an image.
	 * @param workArray byte[] array containing pixel values
	 * @param phase FORE or BACK for foreground of background respectively
	 * @return particleLabels int[] array containing label associating every pixel
	 *         with a particle
	 */
	private int[][] firstIDAttribution(final ImagePlus imp,
		final byte[][] workArray, final int phase)
	{
		final long startTime = System.nanoTime();
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int nSlices = imp.getImageStackSize();
		final int wh = w * h;
		final int nProcessors = Runtime.getRuntime().availableProcessors();
		final int minSlicesPerChunk = 10;
		
		//set up number of chunks
		final int nChunks =
				nSlices < minSlicesPerChunk * nProcessors ?
				(int) Math.ceil((double) nSlices / (double) minSlicesPerChunk) :
				nProcessors;
		
		//set up chunk sizes - last chunk is the remainder
		final int slicesPerChunk = (int) Math.ceil((double) nSlices / (double) nChunks);
		
		//set up start slice array
		final int[] startSlices = new int[nChunks];
		for (int i = 0; i < nChunks; i++) {
			startSlices[i] = i * slicesPerChunk;
		}
		
		
		//set up label offsets to avoid collisions between chunks
		final int chunkLabelSpace = MAX_LABEL / nChunks;
		final int[] chunkIDOffsets = new int[nChunks];
		for (int i = 0; i < nChunks; i++) {
			chunkIDOffsets[i] = i * chunkLabelSpace;
		}
			
		//set up a map split into one per chunk
		final ArrayList<ArrayList<HashSet<Integer>>> chunkMaps = new ArrayList<>(nChunks);
		for (int chunk = 0; chunk < nChunks; chunk++) {
			//assume there is a new particle label for every 10000 pixels
			final int initialArrayCapacity = 1 + w * h * slicesPerChunk / 10000;
			final ArrayList<HashSet<Integer>> map = new ArrayList<>(initialArrayCapacity);
			final int initialHashSetCapacity = 1;
			final int IDoffset = chunkIDOffsets[chunk];
			for (int j = 0; j < initialArrayCapacity; j++) {
				//create a new set containing a single value of j + IDoffset (the root) 
				final HashSet<Integer> set = new HashSet<>(initialHashSetCapacity);
				set.add(j + IDoffset);
				map.add(set);
			}
			chunkMaps.add(map);
		}
		
		//set up the particle label stack
		final int[][] particleLabels = new int[nSlices][wh];
		
		//set up the threads (one thread per chunk)
		final Thread[] threads = new Thread[nChunks];
		
		for (int thread = 0; thread < nChunks; thread++) {
			//each chunk is processed in a single thread
			final int chunk = thread;
			//the starting ID for each chunk is the offset
			final int IDoffset = chunkIDOffsets[chunk]; 
			threads[chunk] = new Thread(() -> {
				//get the Array of HashSets that relate to this image chunk
				final ArrayList<HashSet<Integer>> chunkMap = chunkMaps.get(chunk);
				
				//label image IDs have the chunk ID offset
				int ID = IDoffset;
				
				if (ID == 0) ID = 1;
				
				final int startSlice = startSlices[chunk];
								
				//final slice of the chunk is the next chunk's start slice minus one for all but the last chunk
				final int endSlice = chunk + 1 < nChunks ? startSlices[chunk + 1] - 1 : nSlices - 1;
								
				if (phase == FORE) {
				  //first slice of the chunk - use 4 neighbourhood to not
					//bleed into prior chunk
					final int[] sliceNbh = new int[4];
					for (int y = 0; y < h; y++) {
						final int rowIndex = y * w;
						for (int x = 0; x < w; x++) {
							final int arrayIndex = rowIndex + x;
							if (workArray[startSlice][arrayIndex] == FORE) {
								// Find the minimum particleLabel in the
								// neighbouring pixels
								get4Neighborhood(sliceNbh, particleLabels, x, y, startSlice, w, h, nSlices);
								
								final int minTag = getMinTag(sliceNbh, ID);
								
								//add neighbourhood to map
								addNeighboursToMap(chunkMap, sliceNbh, minTag, IDoffset);
																
								// assign the smallest particle label from the
								// neighbours to the pixel
								particleLabels[startSlice][arrayIndex] = minTag;
								
								// increment the particle label
								if (minTag == ID) {
									ID++;
									expandMap(chunkMap, ID, IDoffset);
								}
							}
						}
					}
					
					//use 13 neighbourhood for all but first slice
					final int[] nbh = new int[13];
					for (int z = startSlice + 1; z <= endSlice ; z++) {
						for (int y = 0; y < h; y++) {
							final int rowIndex = y * w;
							for (int x = 0; x < w; x++) {
								final int arrayIndex = rowIndex + x;
								if (workArray[z][arrayIndex] == FORE) {

									// Find the minimum particleLabel in the
									// neighbouring pixels
									get13Neighborhood(nbh, particleLabels, x, y, z, w, h, nSlices);
									
									final int minTag = getMinTag(nbh, ID);
									
								  //add neighbourhood to map
									addNeighboursToMap(chunkMap, nbh, minTag, IDoffset);
									
									// assign the smallest particle label from the
									// neighbours to the pixel
									particleLabels[z][arrayIndex] = minTag;
									// increment the particle label
									if (minTag == ID) {
										ID++;
										expandMap(chunkMap, ID, IDoffset);
									}
								}
							}
						}
					}
				}
				
				else if (phase == BACK) {
				  //first slice of the chunk - use 2 neighbourhood to not
					//bleed into prior chunk
					final int[] sliceNbh = new int[2];
					for (int y = 0; y < h; y++) {
						final int rowIndex = y * w;
						for (int x = 0; x < w; x++) {
							final int arrayIndex = rowIndex + x;
							if (workArray[startSlice][arrayIndex] == BACK) {
								// Find the minimum particleLabel in the
								// neighbouring pixels
								get2Neighborhood(sliceNbh, particleLabels, x, y, startSlice, w, h, nSlices);
								
								final int minTag = getMinTag(sliceNbh, ID);
								
								//add neighbourhood to map
								addNeighboursToMap(chunkMap, sliceNbh, minTag, IDoffset);
																
								// assign the smallest particle label from the
								// neighbours to the pixel
								particleLabels[startSlice][arrayIndex] = minTag;
								// increment the particle label
								if (minTag == ID) {
									ID++;
									expandMap(chunkMap, ID, IDoffset);
								}
							}
						}
					}

					//use 3-neighbourhood for all but the first slice
					final int[] nbh = new int[3];
					for (int z = 0; z < nSlices; z++) {
						for (int y = 0; y < h; y++) {
							final int rowIndex = y * w;
							for (int x = 0; x < w; x++) {
								final int arrayIndex = rowIndex + x;
								if (workArray[z][arrayIndex] == BACK) {

									// Find the minimum particleLabel in the
									// neighbouring pixels
									get3Neighborhood(nbh, particleLabels, x, y, z, w, h, nSlices);

									final int minTag = getMinTag(nbh, ID);
									
									addNeighboursToMap(chunkMap, nbh, minTag, IDoffset);
									
									// assign the smallest particle label from the
									// neighbours to the pixel
									particleLabels[z][arrayIndex] = minTag;
									// increment the particle label
									if (minTag == ID) {
										ID++;
										expandMap(chunkMap, ID, IDoffset);
									}
								}
							}
						}
					}
				}
//				//there is always one too many IDs per chunk, so trim the last one off
				chunkMap.remove(chunkMap.size() - 1);
			});
		}
		Multithreader.startAndJoin(threads);
		
		//find neighbours in the previous chunk
		//this will result in occasional HashSet values less than 
		//the chunk's IDoffset, which indicate linkage between chunks
		final Thread[] stitchingThreads = new Thread[nChunks];
		for (int thread = 0; thread < nChunks; thread++) {
			final int chunk = thread;
			stitchingThreads[thread] = new Thread(() -> {
				
				//need only one z per thread
				final int z = startSlices[chunk];
				final ArrayList<HashSet<Integer>> chunkMap = chunkMaps.get(chunk);
				final int IDoffset = chunkIDOffsets[chunk];
				
				if (chunk > 0) {
				if (phase == FORE) {
					final int[] nbh = new int[9];
					for (int y = 0; y < h; y++) {
						final int rowIndex = y * w;
						for (int x = 0; x < w; x++) {
							final int arrayIndex = rowIndex + x;
							if (workArray[z][arrayIndex] == FORE) {
								final int label = particleLabels[z][arrayIndex];
								get9Neighborhood(nbh, particleLabels, x, y, z, w, h, nSlices);
								addChunkNeighboursToMap(chunkMap, nbh, label - IDoffset);
							}
						}
					}
				}

				if (phase == BACK) {
					final int[] nbh = new int[1];
					for (int y = 0; y < h; y++) {
						final int rowIndex = y * w;
						for (int x = 0; x < w; x++) {
							final int arrayIndex = rowIndex + x;
							if (workArray[z][arrayIndex] == FORE) {
								final int label = particleLabels[z][arrayIndex];
								get1Neighborhood(nbh, particleLabels, x, y, z, w);
								addChunkNeighboursToMap(chunkMap, nbh, label - IDoffset);
							}
						}
					}
				}
				}
			});
		}
		Multithreader.startAndJoin(stitchingThreads);
				
		final long labellingCompleteTime = System.nanoTime();
		IJ.log("First labelling complete in "+((labellingCompleteTime - startTime)/1000000)+" ms");
		
 	    //snowball the HashSets, handling the chunk offsets and indexes
		//iterate backwards through the chunk maps
		
		boolean somethingChanged = true;
		while (somethingChanged) {
			somethingChanged = false;
			for (int chunk = nChunks - 1; chunk >= 0 ; chunk--) {
				final ArrayList<HashSet<Integer>> map = chunkMaps.get(chunk);
				final int priorChunk = chunk > 0 ? chunk - 1 : 0;
				final ArrayList<HashSet<Integer>> priorMap = chunkMaps.get(priorChunk);
				final int IDoffset = chunkIDOffsets[chunk];
				final int priorIDoffset = chunkIDOffsets[priorChunk];
				for (int i = map.size() - 1; i >= 0; i--) {
					final HashSet<Integer> set = map.get(i);
					if (!set.isEmpty()) {
						//find the minimum label in the set
						int minLabel = Integer.MAX_VALUE;
						for (Integer label : set) {
							if (label < minLabel)
								minLabel = label;
						}
						//if minimum label is less than this chunk's offset, need
						//to move set to previous chunk's map
						if (minLabel < IDoffset) {
							priorMap.get(minLabel - priorIDoffset).addAll(set);
							set.clear();
							somethingChanged = true;
							continue;
						}
						//move whole set's contents to a lower position in the map
						if (minLabel < i + IDoffset) {
							map.get(minLabel - IDoffset).addAll(set);
							set.clear();
							somethingChanged = true;
							continue;
						}
					}
				}
			}
		}
		final long snowballingCompleteTime = System.nanoTime();
		IJ.log("Snowballing complete in "+((snowballingCompleteTime - labellingCompleteTime)/1000000)+" ms");

		//count unique labels and particles
		int labelCount = 0;
		nParticles = 0;
		for (ArrayList<HashSet<Integer>> map : chunkMaps) {
			for (HashSet<Integer> set : map) {
				if (!set.isEmpty()) {
					labelCount += set.size();
					nParticles++;
				}
			}
		}
		
		//set up a 1D HashMap of HashSets with the minimum label
		//set as the 'root' (key) of the hashMap
		HashMap<Integer, HashSet<Integer>> hashMap = new HashMap<>(labelCount);
		for (ArrayList<HashSet<Integer>> map : chunkMaps) {
			for (HashSet<Integer> set : map) {
				int root = Integer.MAX_VALUE;
				for (Integer label : set) {
					if (label < root)
						root = label;
				}
				hashMap.put(root, set);
			}
		}
		
		//set up a LUT to keep track of the minimum replacement value for each label
		final HashMap<Integer, Integer> lutMap = new HashMap<>(labelCount);
		for (ArrayList<HashSet<Integer>> map : chunkMaps) {
			for (HashSet<Integer> set : map) {
				for (Integer label : set)
					//start so that each label looks up itself
					lutMap.put(label, label);
			}
		}

		//check the hashMap for duplicate appearances and merge sets downwards
		somethingChanged = true;
		while (somethingChanged) {
			somethingChanged = false;
			Iterator<Map.Entry<Integer, HashSet<Integer>>> it = hashMap.entrySet().iterator();
			while (it.hasNext()) {
				Map.Entry<Integer, HashSet<Integer>> pair = it.next();
				HashSet<Integer> set = pair.getValue();
				int key = pair.getKey();
				for (Integer label : set) {
					int lutValue = lutMap.get(label);
					//lower the lut lookup value to the root of this set
					if (lutValue > key) {
						lutMap.put(label, key);
						somethingChanged = true;
					}
					//looks like there is a value in the wrong place
					if (lutValue < key) {
						//move all the set's labels to the lower root
						hashMap.get(lutValue).addAll(set);
						//update all the set's lut lookups with the new root
						for (Integer l : set) {
							lutMap.put(l, lutValue);
						}
						set.clear();
						somethingChanged = true;
						break;
					}
				}
			}
		}
		
		//count number of unique values in the LUT
		HashSet<Integer> lutValues = new HashSet<>();
		Iterator<Map.Entry<Integer, Integer>> itL = lutMap.entrySet().iterator();
		while (itL.hasNext()) {
			Map.Entry<Integer, Integer> pair = itL.next();
			lutValues.add(pair.getValue());
		}
		final int nLabels = lutValues.size();
		
		final long hashMappingCompleteTime = System.nanoTime();
		IJ.log("Hashmapping complete in "+((hashMappingCompleteTime - snowballingCompleteTime)/1000000)+" ms");
		
		//assign incremental replacement values
		//translate old 
		final HashMap<Integer, Integer> lutLut = new HashMap<>(nLabels);
		int value = 1;
		for (Integer lutValue : lutValues) {
			if (lutValue == 0) {
				lutLut.put(0, 0);
				continue;
			}
			lutLut.put(lutValue, value);
			value++;
		}
		
		//lutLut now contains mapping from the old lut value (the lutLut 'key') to the
		//new lut value (lutLut 'value')
		
		Iterator<Map.Entry<Integer, Integer>> itR = lutMap.entrySet().iterator();
		while (itR.hasNext()) {
			Map.Entry<Integer, Integer> pair = itR.next();
			Integer oldLutValue = pair.getValue();
			Integer newLutValue = lutLut.get(oldLutValue);
			pair.setValue(newLutValue);
		}
		
		//translate the HashMap LUT to a chunkwise LUT, to be used in combination
		//with the IDoffsets.
		int[][] lut = new int[nChunks][];
		for (int chunk = 0; chunk < nChunks; chunk++) {
			final int nChunkLabels = chunkMaps.get(chunk).size();
			final int IDoffset = chunkIDOffsets[chunk];
			int[] chunkLut = new int[nChunkLabels];
			for (int i = 0; i < nChunkLabels; i++) {
				chunkLut[i] = lutMap.get(i + IDoffset);
			}
			lut[chunk] = chunkLut;
		}
		
		final long lutCreationTime = System.nanoTime();
		IJ.log("LUT creation complete in "+((lutCreationTime - hashMappingCompleteTime)/1000000)+" ms");
		
		//rewrite the pixel values using the LUT
		applyLUT(particleLabels, lut, chunkIDOffsets, startSlices, w, h, nSlices);
		
		final long lutAppliedTime = System.nanoTime();
		IJ.log("LUT applied in "+((lutAppliedTime - lutCreationTime)/1000000)+" ms");
		
		return particleLabels;
	}
	
		/**
		 * Increase the length of the list of label HashSets to
		 * accommodate the full range of IDs
		 * 
		 * @param map
		 * @param ID
		 * @param IDoffset
		 */
		private static void expandMap(final List<HashSet<Integer>> map,
			final int ID, final int IDoffset) {
			while (ID - IDoffset >= map.size()) {
				final HashSet<Integer> set = new HashSet<>();
				set.add(map.size() + IDoffset);
				map.add(set);
			}
		}
	
	/**
	 * Add all the neighbouring labels of a pixel to the map, except 0
	 * (background) and the pixel's own label, which is already in the map.
	 * 
	 * This chunked version of the map stores label IDs ('centre') in the HashSet and 
	 * uses label ID minus per chunk ID offset as the List index.  
	 * 
	 * In this version the non-zero neighbours' labels  are always bigger than the centre, so 
	 * the centre value is added to the neighbours' map indices.
	 *
	 * @param map a map of LUT values.
	 * @param nbh a neighbourhood in the image.
	 * @param centre current pixel's label (with offset)
	 * @param IDoffset chunk's ID offset
	 */
	private static void addNeighboursToMap(final List<HashSet<Integer>> map,
		final int[] nbh, final int centre, final int IDoffset)
	{
		final int l = nbh.length;
		int lastNonZero = -1;
		for (int i = 0; i < l; i++) {
			final int val = nbh[i];

			// skip background, self-similar, and the last label added
			// adding them again is a redundant waste of time
			if (val == 0 || val == centre || val == lastNonZero)
				continue;
			map.get(val - IDoffset).add(centre);
			lastNonZero = val;
		}
	}
	
	/**
	 * Add all the neighbouring labels of a pixel to the map, except 0
	 * (background). The
	 * LUT gets updated with the minimum neighbour found, but this is only within
	 * the first neighbours and not the minimum label in the pixel's neighbour
	 * network
	 *
	 * @param map a map of LUT values.
	 * @param nbh a neighbourhood in the image.
	 * @param centre current pixel's map index (label - IDoffset)
	 */
	private static void addChunkNeighboursToMap(final List<HashSet<Integer>> map,
		final int[] nbh, final int centre)
	{
		final int l = nbh.length;
		final HashSet<Integer> set = map.get(centre);
		int lastNonZero = -1;
		for (int i = 0; i < l; i++) {
			final int val = nbh[i];
			// skip background
			// and the last non-zero value (already added)
			if (val == 0 || val == lastNonZero)
				continue;
			set.add(val);
			lastNonZero = val;
		}
	}
	
	/**
	 * Apply the LUT in multiple threads
	 * 
	 * @param particleLabels
	 * @param lut
	 * @param w
	 * @param h
	 * @param d
	 */
	private static void applyLUT(final int[][] particleLabels,
		final int[][] lut, final int[] chunkIDOffsets, final int[] startSlices,
		final int w, final int h, final int d)
	{
		final int nChunks = chunkIDOffsets.length;
		
		final Thread[] threads = new Thread[nChunks];
		for (int thread = 0; thread < nChunks; thread++) {
			final int chunk = thread;
			threads[thread] = new Thread(() -> {
				final int startSlice = startSlices[chunk];
				final int endSlice = chunk + 1 < nChunks ? startSlices[chunk + 1] - 1 : d - 1;
				final int IDoffset = chunkIDOffsets[chunk];
				final int[] chunkLut = lut[chunk]; 
				for (int z = startSlice; z <= endSlice; z++) {
					final int[] slice = particleLabels[z];
					final int l = slice.length;
					for (int i = 0; i < l; i++) {
						final int label = slice[i];
						if (label == 0) continue;
						slice[i] = chunkLut[label - IDoffset];
					}
				}
			});
		}
		Multithreader.startAndJoin(threads);
	}

	/**
	 * Get 13 neighborhood of a pixel in a 3D image (0 border conditions)
	 * Longhand, hard-coded for speed. This neighbourhood contains the 
	 * set of pixels that have already been visited by the cursor
	 * as it raster scans in an x-y-z order.
	 *
	 * @param neighborhood a neighbourhood in the image.
	 * @param image 3D image (int[][])
	 * @param x x- coordinate
	 * @param y y- coordinate
	 * @param z z- coordinate (in image stacks the indexes start at 1)
	 * @param w width of the image.
	 * @param h height of the image.
	 * @param d depth of the image.
	 */
	private static void get13Neighborhood(final int[] neighborhood,
		final int[][] image, final int x, final int y, final int z, final int w,
		final int h, final int d)
	{
		final int xm1 = x - 1;
		final int xp1 = x + 1;
		final int ym1 = y - 1;
		final int yp1 = y + 1;
		final int zm1 = z - 1;

		neighborhood[0] = getPixel(image, xm1, ym1, zm1, w, h, d);
		neighborhood[1] = getPixel(image, x, ym1, zm1, w, h, d);
		neighborhood[2] = getPixel(image, xp1, ym1, zm1, w, h, d);
		
		neighborhood[3] = getPixel(image, xm1, y, zm1, w, h, d);
		neighborhood[4] = getPixel(image, x, y, zm1, w, h, d);
		neighborhood[5] = getPixel(image, xp1, y, zm1, w, h, d);
		
		neighborhood[6] = getPixel(image, xm1, yp1, zm1, w, h, d);
		neighborhood[7] = getPixel(image, x, yp1, zm1, w, h, d);
		neighborhood[8] = getPixel(image, xp1, yp1, zm1, w, h, d);
		
		neighborhood[9] = getPixel(image, xm1, ym1, z, w, h, d);
		neighborhood[10] = getPixel(image, x, ym1, z, w, h, d);
		neighborhood[11] = getPixel(image, xp1, ym1, z, w, h, d);
		
		neighborhood[12] = getPixel(image, xm1, y, z, w, h, d);
	}
	
	/**
	 * Get 9 neighborhood of a pixel in a 3D image (0 border conditions)
	 * Longhand, hard-coded for speed. This neighbourhood contains the 
	 * set of pixels in previous plane (z-1) of the pixel's 26-neighbourhood
	 *
	 * @param neighborhood a neighbourhood in the image.
	 * @param image 3D image (int[][])
	 * @param x x- coordinate
	 * @param y y- coordinate
	 * @param z z- coordinate (in image stacks the indexes start at 1)
	 * @param w width of the image.
	 * @param h height of the image.
	 * @param d depth of the image.
	 */
	private static void get9Neighborhood(final int[] neighborhood,
		final int[][] image, final int x, final int y, final int z, final int w,
		final int h, final int d)
	{
		final int xm1 = x - 1;
		final int xp1 = x + 1;
		final int ym1 = y - 1;
		final int yp1 = y + 1;
		final int zm1 = z - 1;

		neighborhood[0] = getPixel(image, xm1, ym1, zm1, w, h, d);
		neighborhood[1] = getPixel(image, x, ym1, zm1, w, h, d);
		neighborhood[2] = getPixel(image, xp1, ym1, zm1, w, h, d);
		
		neighborhood[3] = getPixel(image, xm1, y, zm1, w, h, d);
		neighborhood[4] = getPixel(image, x, y, zm1, w, h, d);
		neighborhood[5] = getPixel(image, xp1, y, zm1, w, h, d);
		
		neighborhood[6] = getPixel(image, xm1, yp1, zm1, w, h, d);
		neighborhood[7] = getPixel(image, x, yp1, zm1, w, h, d);
		neighborhood[8] = getPixel(image, xp1, yp1, zm1, w, h, d);
	}
	
	
	/**
	 * Get 4 neighborhood of a pixel in a 3D image (0 border conditions)
	 * Longhand, hard-coded for speed. This neighbourhood contains the 
	 * set of pixels that have already been visited by the cursor
	 * in the current plane as it raster scans in an x-y order.
	 *
	 * @param neighborhood a neighbourhood in the image.
	 * @param image 3D image (int[][])
	 * @param x x- coordinate
	 * @param y y- coordinate
	 * @param z z- coordinate (in image stacks the indexes start at 1)
	 * @param w width of the image.
	 */
	private static void get4Neighborhood(final int[] neighborhood,
		final int[][] image, final int x, final int y, final int z, final int w, final int h, final int d)
	{
		final int xm1 = x - 1;
		final int xp1 = x + 1;
		final int ym1 = y - 1;
				
		neighborhood[0] = getPixel(image, xm1, ym1, z, w, h, d);
		neighborhood[1] = getPixel(image, x, ym1, z, w, h, d);
		neighborhood[2] = getPixel(image, xp1, ym1, z, w, h, d);
		
		neighborhood[3] = getPixel(image, xm1, y, z, w, h, d);
	}
	
	private static void get3Neighborhood(final int[] neighborhood,
		final int[][] image, final int x, final int y, final int z, final int w,
		final int h, final int d)
	{
		neighborhood[0] = getPixel(image, x - 1, y, z, w, h, d);
		neighborhood[1] = getPixel(image, x, y - 1, z, w, h, d);
		neighborhood[2] = getPixel(image, x, y, z - 1, w, h, d);
	}
	
	private static void get2Neighborhood(final int[] neighborhood,
		final int[][] image, final int x, final int y, final int z, final int w,
		final int h, final int d)
	{
		neighborhood[0] = getPixel(image, x - 1, y, z, w, h, d);
		neighborhood[1] = getPixel(image, x, y - 1, z, w, h, d);
	}
	
	private static void get1Neighborhood(final int[] neighborhood,
		final int[][] image, final int x, final int y, final int z, final int w)
	{
		neighborhood[0] = image[z - 1][x + y * w];
	}


	/**
	 * Get pixel in 3D image (0 border conditions)
	 *
	 * @param image 3D image
	 * @param x x- coordinate
	 * @param y y- coordinate
	 * @param z z- coordinate (in image stacks the indexes start at 1)
	 * @param w width of the image.
	 * @param h height of the image.
	 * @param d depth of the image.
	 * @return corresponding pixel (0 if out of image)
	 */
	private static int getPixel(final int[][] image, final int x, final int y,
		final int z, final int w, final int h, final int d)
	{
		if (withinBounds(x, y, z, w, h, d)) return image[z][x + y * w];

		return 0;
	}
	
	/**
	 * checks whether a pixel at (m, n, o) is within the image boundaries
	 * 
	 * 26- and 6-neighbourhood version 
	 * 
	 * @param m x coordinate
	 * @param n y coordinate
	 * @param o z coordinate
	 * @param w image width
	 * @param h image height
	 * @param d image depth
	 * @return true if the pixel is within the image bounds
	 */
	private static boolean withinBounds(final int m, final int n, final int o,
		final int w, final int h, final int d)
	{
		return (m >= 0 && m < w && n >= 0 && n < h && o >= 0 && o < d);
	}	

	private static int getMinTag(final int[] neighbourhood, final int ID) {
		final int l = neighbourhood.length;
		int minTag = ID;
		for (int i = 0; i < l; i++) {
			final int tagv = neighbourhood[i];
			if (tagv == 0) continue;
			if (tagv < minTag) minTag = tagv;
		}
		return minTag;
	}
}