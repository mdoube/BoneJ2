package org.bonej.plugins;

import java.util.concurrent.atomic.AtomicInteger;

import org.bonej.util.Multithreader;

import ij.ImagePlus;
import ij.measure.Calibration;

public class ParticleAnalysis {

	public ParticleAnalysis() {

	}

	/**
	 * Remove edge-touching, too small and too big particles.
	 * 
	 * 
	 * Relabels the particles to be continuous from 1 to maxParticle and updates
	 * particleSizes accordingly. workArray is updated too.
	 * 
	 * @param imp            Input image. Needed for calibration
	 * @param particleLabels Particle label image array
	 * @param workArray      Binary work array
	 * @param particleSizes  List of particle sizes in pixels indexed by particle ID
	 * @param phase          Foreground or background
	 * @param doExclude      true to remove all particles touching a side
	 * @param min            minimum volume in calibrated units to include
	 * @param max            minimum volume in calibrated units to include
	 */
	public void filterParticles(final ImagePlus imp, final int[][] particleLabels, final byte[][] workArray,
			long[] particleSizes, final int phase, final boolean doExclude, final double min, final double max) {

		// flag to check whether sizes & labels arrays need to be updated
		boolean runLutNeeded = false;

		int nParticles = particleSizes.length;

		// set up the replacement lut
		final int[] lut = new int[nParticles];
		for (int i = 0; i < nParticles; i++)
			lut[i] = i;

		if (min > 0 || max < Double.POSITIVE_INFINITY) {
			// do the size filtering check
			Calibration cal = imp.getCalibration();
			final double pxVol = cal.pixelDepth * cal.pixelHeight * cal.pixelWidth;

			for (int i = 0; i < nParticles; i++) {
				final double particleVolume = particleSizes[i] * pxVol;
				// if this particle is outside min & max limits
				if (particleVolume > max || particleVolume < min) {
					// set lut value to 0 schedules for deletion
					lut[i] = 0;
					runLutNeeded = true;
				}
			}
		}

		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int d = imp.getImageStackSize();

		if (doExclude) {
			// do the edge filtering check. Set all edge particles to 0 in the lut

			// scan faces
			// top and bottom faces
			for (int y = 0; y < h; y++) {
				final int index = y * w;
				for (int x = 0; x < w; x++) {
					final int pt = particleLabels[0][index + x];
					if (pt > 0) {
						lut[pt] = 0;
						runLutNeeded = true;
					}
					final int pb = particleLabels[d - 1][index + x];
					if (pb > 0) {
						lut[pb] = 0;
						runLutNeeded = true;
					}
				}
			}

			// west and east faces
			for (int z = 0; z < d; z++) {
				for (int y = 0; y < h; y++) {
					final int pw = particleLabels[z][y * w];
					final int pe = particleLabels[z][y * w + w - 1];
					if (pw > 0) {
						lut[pw] = 0;
						runLutNeeded = true;
					}
					if (pe > 0) {
						lut[pe] = 0;
						runLutNeeded = true;
					}
				}
			}

			// north and south faces
			final int lastRow = w * (h - 1);
			for (int z = 0; z < d; z++) {
				for (int x = 0; x < w; x++) {
					final int pn = particleLabels[z][x];
					final int ps = particleLabels[z][lastRow + x];
					if (pn > 0) {
						lut[pn] = 0;
						runLutNeeded = true;
					}
					if (ps > 0) {
						lut[ps] = 0;
						runLutNeeded = true;
					}
				}
			}
		}

		// check the arrays only if needed
		if (runLutNeeded) {

			// minimise the lut and count non-zeros
			int nonZeroCount = 0;
			for (int i = 0; i < nParticles; i++) {
				final int lutValue = lut[i];
				if (lutValue != 0) {
					nonZeroCount++;
					lut[i] = nonZeroCount;
				}
			}
			// lut is now 0 for particles to be deleted or a lower value to get rid of gaps

			// reset nParticles, +1 is for particle 0 (background particle)
			nParticles = nonZeroCount + 1;

			// replace particle sizes based on lut
			long[] filteredParticleSizes = new long[nParticles];

			final int l = particleSizes.length;
			for (int i = 0; i < l; i++) {
				final long size = particleSizes[i];
				final int lutValue = lut[i];
				if (lutValue != 0) {
					filteredParticleSizes[lutValue] = size;
				}
			}

			// particleSizes now has the shorter length to match the new particleLabels
			particleSizes = filteredParticleSizes;

			// replace labels based on lut

			// handle both phases in the workArray
			final byte flip;
			if (phase == ConnectedComponents.FORE) {
				flip = 0;
			} else {
				flip = (byte) 255;
			}

			final int wh = w * h;

			AtomicInteger ai = new AtomicInteger(0);

			final Thread[] threads = Multithreader.newThreads();
			for (int thread = 0; thread < threads.length; thread++) {
				threads[thread] = new Thread(() -> {
					for (int z = ai.getAndIncrement(); z < d; z = ai.getAndIncrement()) {
						final int[] particleLabelSlice = particleLabels[z];
						final byte[] workArraySlice = workArray[z];
						for (int i = 0; i < wh; i++) {
							final int oldLabel = particleLabelSlice[i];
							if (oldLabel == 0)
								continue;
							final int newLabel = lut[oldLabel];
							particleLabelSlice[i] = newLabel;
							if (newLabel == 0)
								workArraySlice[i] = flip;
						}
					}
				});
			}
			Multithreader.startAndJoin(threads);
		}
	}
}
