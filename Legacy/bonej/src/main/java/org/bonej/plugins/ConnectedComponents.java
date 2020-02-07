package org.bonej.plugins;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import org.bonej.util.Multithreader;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

public class ConnectedComponents {

	/** Foreground value */
	static final int FORE = -1;
	/** Background value */
	static final int BACK = 0;
	/** 2^23 - greatest integer that can be represented precisely by a float */
	static final int MAX_LABEL = 8388608;

	/** number of particle labels */
	private static int nParticles;

	private static byte[][] workArray;

	public ConnectedComponents() {

	}

	/**
	 * Run connected components filter on a binary image
	 * 
	 * @param imp   Input ImagePlus, must be 2D or 3D and binary (0 & 255)
	 * @param phase either foreground or background
	 * @return 2D int array with the same dimensions as the input image, with
	 *         individual connected components labelled with a unique, consecutive
	 *         label.
	 */
	int[][] run(final ImagePlus imp, final int phase) {
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final int nSlices = imp.getImageStackSize();
		final int nProcessors = Runtime.getRuntime().availableProcessors();
		final int minSlicesPerChunk = 10;

		// set up number of chunks
		final int nChunks = nSlices < minSlicesPerChunk * nProcessors
				? (int) Math.ceil((double) nSlices / (double) minSlicesPerChunk)
				: nProcessors;

		// set up chunk sizes - last chunk is the remainder
		final int slicesPerChunk = (int) Math.ceil((double) nSlices / (double) nChunks);

		// set up start slice array
		final int[] startSlices = new int[nChunks];
		for (int i = 0; i < nChunks; i++) {
			startSlices[i] = i * slicesPerChunk;
		}

		// set up label offsets to avoid collisions between chunks
		final int chunkLabelSpace = MAX_LABEL / nChunks;
		final int[] chunkIDOffsets = new int[nChunks];
		for (int i = 0; i < nChunks; i++) {
			chunkIDOffsets[i] = i * chunkLabelSpace;
		}

		// set up a map split into one per chunk
		final ArrayList<ArrayList<HashSet<Integer>>> chunkMaps = new ArrayList<>(nChunks);
		for (int chunk = 0; chunk < nChunks; chunk++) {
			// assume there is a new particle label for every 10000 pixels
			final int initialArrayCapacity = 1 + w * h * slicesPerChunk / 10000;
			final ArrayList<HashSet<Integer>> map = new ArrayList<>(initialArrayCapacity);
			chunkMaps.add(map);
		}

		// set up the work array
		makeWorkArray(imp);

		//do a first labelling and map first degree neighbours
		int[][] particleLabels = firstIDAttribution(chunkMaps, chunkIDOffsets, startSlices, w, h, nSlices, phase);

		//merge neighbour networks and generate a LUT
		final int[][] lut = generateLut(chunkMaps, chunkIDOffsets);
		
		// rewrite the pixel values using the LUT
		applyLUT(particleLabels, lut, chunkIDOffsets, startSlices, w, h, nSlices);

		return particleLabels;
	}

	private static int[][] generateLut(ArrayList<ArrayList<HashSet<Integer>>> chunkMaps, int[] chunkIDOffsets) {
		// snowball the HashSets, handling the chunk offsets and indexes
		snowball(chunkMaps, chunkIDOffsets);

		HashMap<Integer, Integer> lutMap = makeLutMap(chunkMaps);

		return lutFromLutMap(lutMap, chunkMaps, chunkIDOffsets);
	}

	/**
	 * Create a work array and store it as a field of this instance, which can be
	 * retrieved with getWorkArray
	 *
	 * @param imp an image.
	 */
	static void makeWorkArray(final ImagePlus imp) {
		final int s = imp.getStackSize();
		final int p = imp.getWidth() * imp.getHeight();
		workArray = new byte[s][p];
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
	}

	/**
	 * Go through all pixels and assign initial particle label.
	 *
	 * @param chunkMaps an image.
	 * @param phase     FORE or BACK for foreground of background respectively
	 * @return particleLabels int[] array containing label associating every pixel
	 *         with a particle
	 */
	private static int[][] firstIDAttribution(final ArrayList<ArrayList<HashSet<Integer>>> chunkMaps,
			final int[] chunkIDOffsets, final int[] startSlices, final int w, final int h, final int nSlices,
			final int phase) {

		final int nChunks = chunkIDOffsets.length;
		final int wh = w * h;
		// set up the particle label stack
		final int[][] particleLabels = new int[nSlices][wh];

		// set up the threads (one thread per chunk)
		final Thread[] threads = new Thread[nChunks];

		for (int thread = 0; thread < nChunks; thread++) {
			// each chunk is processed in a single thread
			final int chunk = thread;
			// the starting ID for each chunk is the offset
			final int IDoffset = chunkIDOffsets[chunk];
			threads[chunk] = new Thread(() -> {
				// get the Array of HashSets that relate to this image chunk
				final ArrayList<HashSet<Integer>> chunkMap = chunkMaps.get(chunk);

				// label image IDs have the chunk ID offset
				int ID = IDoffset;

				if (ID == 0)
					ID = 1;

				final int startSlice = startSlices[chunk];

				// final slice of the chunk is the next chunk's start slice minus one for all
				// but the last chunk
				final int endSlice = chunk + 1 < nChunks ? startSlices[chunk + 1] - 1 : nSlices - 1;

				if (phase == FORE) {
					// first slice of the chunk - use 4 neighbourhood to not
					// bleed into prior chunk
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

								// add neighbourhood to map
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

					// use 13 neighbourhood for all but first slice
					final int[] nbh = new int[13];
					for (int z = startSlice + 1; z <= endSlice; z++) {
						for (int y = 0; y < h; y++) {
							final int rowIndex = y * w;
							for (int x = 0; x < w; x++) {
								final int arrayIndex = rowIndex + x;
								if (workArray[z][arrayIndex] == FORE) {

									// Find the minimum particleLabel in the
									// neighbouring pixels
									get13Neighborhood(nbh, particleLabels, x, y, z, w, h, nSlices);

									final int minTag = getMinTag(nbh, ID);

									// add neighbourhood to map
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
					// first slice of the chunk - use 2 neighbourhood to not
					// bleed into prior chunk
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

								// add neighbourhood to map
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

					// use 3-neighbourhood for all but the first slice
					final int[] nbh = new int[3];
					for (int z = startSlice + 1; z <= endSlice; z++) {
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
				// there is always one too many IDs per chunk, so trim the last one off
				chunkMap.remove(chunkMap.size() - 1);
			});
		}
		Multithreader.startAndJoin(threads);

		// find neighbours in the previous chunk
		// this will result in occasional HashSet values less than
		// the chunk's IDoffset, which indicate linkage between chunks
		final Thread[] stitchingThreads = new Thread[nChunks];
		for (int thread = 0; thread < nChunks; thread++) {
			final int chunk = thread;
			stitchingThreads[thread] = new Thread(() -> {

				// need only one z per thread
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
								if (workArray[z][arrayIndex] == BACK) {
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

		return particleLabels;
	}

	/**
	 * 
	 * @param chunkMaps
	 * @param chunkIDOffsets
	 * @param nChunks
	 */
	private static void snowball(final ArrayList<ArrayList<HashSet<Integer>>> chunkMaps, final int[] chunkIDOffsets) {
		// iterate backwards through the chunk maps

		final int nChunks = chunkIDOffsets.length;

		for (int chunk = nChunks - 1; chunk >= 0; chunk--) {
			final ArrayList<HashSet<Integer>> map = chunkMaps.get(chunk);
			final int priorChunk = chunk > 0 ? chunk - 1 : 0;
			final ArrayList<HashSet<Integer>> priorMap = chunkMaps.get(priorChunk);
			final int IDoffset = chunkIDOffsets[chunk];
			final int priorIDoffset = chunkIDOffsets[priorChunk];
			for (int i = map.size() - 1; i >= 0; i--) {
				final HashSet<Integer> set = map.get(i);
				if (!set.isEmpty()) {
					// find the minimum label in the set
					int minLabel = Integer.MAX_VALUE;
					for (Integer label : set) {
						if (label < minLabel)
							minLabel = label;
					}
					// if minimum label is less than this chunk's offset, need
					// to move set to previous chunk's map
					if (minLabel < IDoffset) {
						priorMap.get(minLabel - priorIDoffset).addAll(set);
						set.clear();
						continue;
					}
					// move whole set's contents to a lower position in the map
					if (minLabel < i + IDoffset) {
						map.get(minLabel - IDoffset).addAll(set);
						set.clear();
						continue;
					}
				}
			}
		}
	}

	/**
	 * 
	 * @param chunkMaps
	 * @return
	 */
	private static HashMap<Integer, Integer> makeLutMap(final ArrayList<ArrayList<HashSet<Integer>>> chunkMaps) {
		// count unique labels and particles
		int labelCount = 0;
		for (ArrayList<HashSet<Integer>> map : chunkMaps) {
			for (HashSet<Integer> set : map) {
				if (!set.isEmpty())
					labelCount += set.size();
			}
		}

		// set up a 1D HashMap of HashSets with the minimum label
		// set as the 'root' (key) of the hashMap
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

		// set up a LUT to keep track of the minimum replacement value for each label
		final HashMap<Integer, Integer> lutMap = new HashMap<>(labelCount);
		for (ArrayList<HashSet<Integer>> map : chunkMaps) {
			for (HashSet<Integer> set : map) {
				for (Integer label : set) {
					// start so that each label looks up itself
					lutMap.put(label, label);
				}
			}
		}

		// check the hashMap for duplicate appearances and merge sets downwards
		boolean somethingChanged = true;
		while (somethingChanged) {
			IJ.log("checking hashMap");
			somethingChanged = false;
			Iterator<Map.Entry<Integer, HashSet<Integer>>> it = hashMap.entrySet().iterator();
			while (it.hasNext()) {
				Map.Entry<Integer, HashSet<Integer>> pair = it.next();
				HashSet<Integer> set = pair.getValue();
				int key = pair.getKey();
				for (Integer label : set) {
					int lutValue = lutMap.get(label);
					// lower the lut lookup value to the root of this set
					if (lutValue > key) {
						lutMap.put(label, key);
						somethingChanged = true;
					}
					// looks like there is a value in the wrong place
					if (lutValue < key) {
						// move all the set's labels to the lower root
						hashMap.get(lutValue).addAll(set);
						// update all the set's lut lookups with the new root
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
		return lutMap;
	}

	/**
	 * 
	 * @param lutMap
	 * @param chunkMaps
	 * @param chunkIDOffsets
	 * @return
	 */
	private static int[][] lutFromLutMap(final HashMap<Integer, Integer> lutMap,
			final ArrayList<ArrayList<HashSet<Integer>>> chunkMaps, final int[] chunkIDOffsets) {
		// count number of unique labels in the LUT
		HashSet<Integer> lutLabels = new HashSet<>();
		Iterator<Map.Entry<Integer, Integer>> itL = lutMap.entrySet().iterator();
		while (itL.hasNext()) {
			Map.Entry<Integer, Integer> pair = itL.next();
			lutLabels.add(pair.getValue());
		}
		final int nLabels = lutLabels.size();
		nParticles = nLabels;

		// assign incremental replacement values
		// translate old
		final HashMap<Integer, Integer> lutLut = new HashMap<>(nLabels);
		int value = 1;
		for (Integer lutValue : lutLabels) {
			if (lutValue == 0) {
				lutLut.put(0, 0);
				continue;
			}
			lutLut.put(lutValue, value);
			value++;
		}

		// lutLut now contains mapping from the old lut value (the lutLut 'key') to the
		// new lut value (lutLut 'value')

		Iterator<Map.Entry<Integer, Integer>> itR = lutMap.entrySet().iterator();
		while (itR.hasNext()) {
			Map.Entry<Integer, Integer> pair = itR.next();
			Integer oldLutValue = pair.getValue();
			Integer newLutValue = lutLut.get(oldLutValue);
			pair.setValue(newLutValue);
		}

		// translate the HashMap LUT to a chunkwise LUT, to be used in combination
		// with the IDoffsets.
		final int nChunks = chunkIDOffsets.length;
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
		return lut;
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
	private static void applyLUT(final int[][] particleLabels, final int[][] lut, final int[] chunkIDOffsets,
			final int[] startSlices, final int w, final int h, final int d) {
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
						if (label == 0)
							continue;
						slice[i] = chunkLut[label - IDoffset];
					}
				}
			});
		}
		Multithreader.startAndJoin(threads);
	}

	/**
	 * Add all the neighbouring labels of a pixel to the map, except 0 (background)
	 * and the pixel's own label, which is already in the map.
	 * 
	 * This chunked version of the map stores label IDs ('centre') in the HashSet
	 * and uses label ID minus per chunk ID offset as the List index.
	 * 
	 * In this version the non-zero neighbours' labels are always bigger than the
	 * centre, so the centre value is added to the neighbours' map indices.
	 *
	 * @param map      a map of LUT values.
	 * @param nbh      a neighbourhood in the image.
	 * @param centre   current pixel's label (with offset)
	 * @param IDoffset chunk's ID offset
	 */
	private static void addNeighboursToMap(final List<HashSet<Integer>> map, final int[] nbh, final int centre,
			final int IDoffset) {
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
	 * Add all the neighbouring labels of a pixel to the map, except 0 (background).
	 * The LUT gets updated with the minimum neighbour found, but this is only
	 * within the first neighbours and not the minimum label in the pixel's
	 * neighbour network
	 *
	 * @param map    a map of LUT values.
	 * @param nbh    a neighbourhood in the image.
	 * @param centre current pixel's map index (label - IDoffset)
	 */
	private static void addChunkNeighboursToMap(final List<HashSet<Integer>> map, final int[] nbh, final int centre) {
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
	 * Get 13 neighborhood of a pixel in a 3D image (0 border conditions) Longhand,
	 * hard-coded for speed. This neighbourhood contains the set of pixels that have
	 * already been visited by the cursor as it raster scans in an x-y-z order.
	 *
	 * @param neighborhood a neighbourhood in the image.
	 * @param image        3D image (int[][])
	 * @param x            x- coordinate
	 * @param y            y- coordinate
	 * @param z            z- coordinate (in image stacks the indexes start at 1)
	 * @param w            width of the image.
	 * @param h            height of the image.
	 * @param d            depth of the image.
	 */
	private static void get13Neighborhood(final int[] neighborhood, final int[][] image, final int x, final int y,
			final int z, final int w, final int h, final int d) {
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
	 * Get 9 neighborhood of a pixel in a 3D image (0 border conditions) Longhand,
	 * hard-coded for speed. This neighbourhood contains the set of pixels in
	 * previous plane (z-1) of the pixel's 26-neighbourhood
	 *
	 * @param neighborhood a neighbourhood in the image.
	 * @param image        3D image (int[][])
	 * @param x            x- coordinate
	 * @param y            y- coordinate
	 * @param z            z- coordinate (in image stacks the indexes start at 1)
	 * @param w            width of the image.
	 * @param h            height of the image.
	 * @param d            depth of the image.
	 */
	private static void get9Neighborhood(final int[] neighborhood, final int[][] image, final int x, final int y,
			final int z, final int w, final int h, final int d) {
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
	 * Get 4 neighborhood of a pixel in a 3D image (0 border conditions) Longhand,
	 * hard-coded for speed. This neighbourhood contains the set of pixels that have
	 * already been visited by the cursor in the current plane as it raster scans in
	 * an x-y order.
	 *
	 * @param neighborhood a neighbourhood in the image.
	 * @param image        3D image (int[][])
	 * @param x            x- coordinate
	 * @param y            y- coordinate
	 * @param z            z- coordinate (in image stacks the indexes start at 1)
	 * @param w            width of the image.
	 */
	private static void get4Neighborhood(final int[] neighborhood, final int[][] image, final int x, final int y,
			final int z, final int w, final int h, final int d) {
		final int xm1 = x - 1;
		final int xp1 = x + 1;
		final int ym1 = y - 1;

		neighborhood[0] = getPixel(image, xm1, ym1, z, w, h, d);
		neighborhood[1] = getPixel(image, x, ym1, z, w, h, d);
		neighborhood[2] = getPixel(image, xp1, ym1, z, w, h, d);

		neighborhood[3] = getPixel(image, xm1, y, z, w, h, d);
	}

	private static void get3Neighborhood(final int[] neighborhood, final int[][] image, final int x, final int y,
			final int z, final int w, final int h, final int d) {
		neighborhood[0] = getPixel(image, x - 1, y, z, w, h, d);
		neighborhood[1] = getPixel(image, x, y - 1, z, w, h, d);
		neighborhood[2] = getPixel(image, x, y, z - 1, w, h, d);
	}

	private static void get2Neighborhood(final int[] neighborhood, final int[][] image, final int x, final int y,
			final int z, final int w, final int h, final int d) {
		neighborhood[0] = getPixel(image, x - 1, y, z, w, h, d);
		neighborhood[1] = getPixel(image, x, y - 1, z, w, h, d);
	}

	private static void get1Neighborhood(final int[] neighborhood, final int[][] image, final int x, final int y,
			final int z, final int w) {
		neighborhood[0] = image[z - 1][x + y * w];
	}

	/**
	 * Get pixel in 3D image (0 border conditions)
	 *
	 * @param image 3D image
	 * @param x     x- coordinate
	 * @param y     y- coordinate
	 * @param z     z- coordinate (in image stacks the indexes start at 1)
	 * @param w     width of the image.
	 * @param h     height of the image.
	 * @param d     depth of the image.
	 * @return corresponding pixel (0 if out of image)
	 */
	private static int getPixel(final int[][] image, final int x, final int y, final int z, final int w, final int h,
			final int d) {
		if (withinBounds(x, y, z, w, h, d))
			return image[z][x + y * w];

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
	private static boolean withinBounds(final int m, final int n, final int o, final int w, final int h, final int d) {
		return (m >= 0 && m < w && n >= 0 && n < h && o >= 0 && o < d);
	}

	private static int getMinTag(final int[] neighbourhood, final int ID) {
		final int l = neighbourhood.length;
		int minTag = ID;
		for (int i = 0; i < l; i++) {
			final int tagv = neighbourhood[i];
			if (tagv == 0)
				continue;
			if (tagv < minTag)
				minTag = tagv;
		}
		return minTag;
	}

	/**
	 * Increase the length of the list of label HashSets to accommodate the full
	 * range of IDs
	 * 
	 * @param map
	 * @param ID
	 * @param IDoffset
	 */
	private static void expandMap(final List<HashSet<Integer>> map, final int ID, final int IDoffset) {
		while (ID - IDoffset >= map.size()) {
			final HashSet<Integer> set = new HashSet<>();
			set.add(map.size() + IDoffset);
			map.add(set);
		}
	}

	public int getNParticles() {
		final int n = nParticles;
		return n;
	}

	public byte[][] getWorkArray() {
		return workArray;
	}

}