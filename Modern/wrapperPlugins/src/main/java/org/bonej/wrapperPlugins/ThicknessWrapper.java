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

package org.bonej.wrapperPlugins;

import static org.bonej.utilities.ImagePlusUtil.cleanDuplicate;
import static org.bonej.wrapperPlugins.CommonMessages.HAS_CHANNEL_DIMENSIONS;
import static org.bonej.wrapperPlugins.CommonMessages.HAS_TIME_DIMENSIONS;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_3D_IMAGE;
import static org.bonej.wrapperPlugins.CommonMessages.NOT_8_BIT_BINARY_IMAGE;
import static org.bonej.wrapperPlugins.CommonMessages.NO_IMAGE_OPEN;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import net.imagej.patcher.LegacyInjector;
import net.imagej.table.DefaultColumn;
import net.imagej.table.Table;

import org.bonej.utilities.ImagePlusUtil;
import org.bonej.utilities.RoiManagerUtil;
import org.bonej.utilities.SharedTable;
import org.bonej.wrapperPlugins.wrapperUtils.Common;
import org.bonej.wrapperPlugins.wrapperUtils.ResultUtils;
import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.ContextCommand;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;
import org.scijava.widget.ChoiceWidget;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.LutLoader;
import ij.plugin.frame.RoiManager;
import ij.process.LUT;
import ij.process.StackStatistics;
import sc.fiji.localThickness.LocalThicknessWrapper;

/**
 * An ImageJ2 command that wraps the sc.fiji.localThickness plugin
 *
 * @author Richard Domander
 */
@Plugin(type = Command.class, menuPath = "Plugins>BoneJ>Thickness")
public class ThicknessWrapper extends ContextCommand {

	static {
		LegacyInjector.preinit();
	}

	@Parameter(validater = "validateImage")
	private ImagePlus inputImage;

	@Parameter(label = "Calculate:",
		description = "Which thickness measures to calculate",
		style = ChoiceWidget.RADIO_BUTTON_VERTICAL_STYLE, choices = {
			"Trabecular thickness", "Trabecular spacing", "Both" })
	private String mapChoice = "Trabecular thickness";

	@Parameter(label = "Show thickness maps",
		description = "Show resulting map images after calculations",
		required = false)
	private boolean showMaps = true;

	@Parameter(label = "Mask thickness maps",
		description = "Remove pixel artifacts from the thickness maps",
		required = false)
	private boolean maskArtefacts = true;

	@Parameter(label = "Crop to ROI manager",
		description = "Limit the maps to the ROIs in the ROI manager",
		persist = false, required = false)
	private boolean cropToRois;

	@Parameter(label = "Trabecular thickness", type = ItemIO.OUTPUT)
	private ImagePlus trabecularMap;

	@Parameter(label = "Trabecular spacing", type = ItemIO.OUTPUT)
	private ImagePlus spacingMap;

	/**
	 * The calculated thickness statistics in a {@link Table}
	 * <p>
	 * Null if there are no results
	 * </p>
	 */
	@Parameter(type = ItemIO.OUTPUT, label = "BoneJ results")
	private Table<DefaultColumn<Double>, Double> resultsTable;

	@Parameter
	private UIService uiService;

	@Parameter
	private StatusService statusService;

	private boolean foreground;
	private LocalThicknessWrapper localThickness;
	private boolean anisotropyWarned;

	@Override
	public void run() {
		final List<Boolean> mapOptions = getMapOptions();
		createLocalThickness();
		final Map<Boolean, ImagePlus> thicknessMaps = new HashMap<>();
		mapOptions.forEach(foreground -> {
			prepareRun(foreground);
			statusService.showStatus("Thickness: creating thickness map");
			final ImagePlus map = createMap();
			statusService.showStatus("Thickness: calculating results");
			addMapResults(map);
			thicknessMaps.put(foreground, map);
		});
		if (SharedTable.hasData()) {
			resultsTable = SharedTable.getTable();
		}
		if (showMaps) {
			LUT fire = Common.makeFire();
			trabecularMap = thicknessMaps.get(true);
			StackStatistics trabecularStats = new StackStatistics(trabecularMap);
			trabecularMap.setDisplayRange(0.0, trabecularStats.max);
			trabecularMap.setLut(fire);
			spacingMap = thicknessMaps.get(false);
			StackStatistics spacingStats = new StackStatistics(spacingMap);
			spacingMap.setDisplayRange(0.0, spacingStats.max);
			spacingMap.setLut(fire);
		}
	}

	private void addMapResults(final ImagePlus map) {
		if (map == null) {
			return;
		}
		final String label = inputImage.getTitle();
		final String unitHeader = ResultUtils.getUnitHeader(map);
		final String prefix = foreground ? "Tb.Th" : "Tb.Sp";
		final StackStatistics resultStats = new StackStatistics(map);
		double mean = resultStats.mean;
		double stdDev = resultStats.stdDev;
		double max = resultStats.max;

		if (resultStats.pixelCount == 0) {
			// All pixels are background (NaN), stats not applicable
			mean = Double.NaN;
			stdDev = Double.NaN;
			max = Double.NaN;
		}

		SharedTable.add(label, prefix + " Mean " + unitHeader, mean);
		SharedTable.add(label, prefix + " Std Dev " + unitHeader, stdDev);
		SharedTable.add(label, prefix + " Max " + unitHeader, max);
	}

	private void createLocalThickness() {
		localThickness = new LocalThicknessWrapper();
		localThickness.setSilence(true);
		localThickness.setShowOptions(false);
		localThickness.maskThicknessMap = maskArtefacts;
		localThickness.calibratePixels = true;
	}

	private ImagePlus createMap() {
		final ImagePlus image;

		if (cropToRois) {
			final RoiManager roiManager = RoiManager.getInstance2();
			final Optional<ImageStack> stackOptional = RoiManagerUtil.cropToRois(
				roiManager, inputImage.getStack(), true, 0x00);
			if (!stackOptional.isPresent()) {
				cancel("Can't crop without valid ROIs in the ROIManager");
				return null;
			}
			image = new ImagePlus(inputImage.getTitle(), stackOptional.get());
		}
		else {
			image = cleanDuplicate(inputImage);
		}

		return localThickness.processImage(image);
	}

	// region -- Helper methods --
	private List<Boolean> getMapOptions() {
		final List<Boolean> mapOptions = new ArrayList<>();
		if ("Trabecular thickness".equals(mapChoice)) {
			mapOptions.add(true);
		}
		else if ("Trabecular spacing".equals(mapChoice)) {
			mapOptions.add(false);
		}
		else if ("Both".equals(mapChoice)) {
			mapOptions.add(true);
			mapOptions.add(false);
		}
		else {
			throw new IllegalArgumentException("Unexpected map choice");
		}
		return mapOptions;
	}

	private void prepareRun(final boolean foreground) {
		this.foreground = foreground;
		final String suffix = foreground ? "_Tb.Th" : "_Tb.Sp";
		localThickness.setTitleSuffix(suffix);
		localThickness.inverse = !foreground;
	}

	@SuppressWarnings("unused")
	private void validateImage() {
		if (inputImage == null) {
			cancel(NO_IMAGE_OPEN);
			return;
		}

		if (!ImagePlusUtil.is3D(inputImage)) {
			cancel(NOT_3D_IMAGE);
			return;
		}

		if (inputImage.getNChannels() > 1) {
			cancel(HAS_CHANNEL_DIMENSIONS + ". Please split the channels.");
			return;
		}
		if (inputImage.getNFrames() > 1) {
			cancel(HAS_TIME_DIMENSIONS + ". Please split the hyperstack.");
			return;
		}

		if (!ImagePlusUtil.isBinaryColour(inputImage) || inputImage
			.getBitDepth() != 8)
		{
			cancel(NOT_8_BIT_BINARY_IMAGE);
			return;
		}

		if (!anisotropyWarned) {
			if (!Common.warnAnisotropy(inputImage, uiService)) {
				cancel(null);
			}
			anisotropyWarned = true;
		}
	}

	// endregion
}
