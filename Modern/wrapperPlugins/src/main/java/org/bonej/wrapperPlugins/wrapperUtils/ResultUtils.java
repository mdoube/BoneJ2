
package org.bonej.wrapperPlugins.wrapperUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.TypedAxis;
import net.imagej.space.AnnotatedSpace;
import net.imagej.table.GenericColumn;
import net.imagej.table.LongColumn;
import net.imagej.units.UnitService;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import org.bonej.utilities.AxisUtils;
import org.bonej.wrapperPlugins.wrapperUtils.HyperstackUtils.Subspace;
import org.scijava.util.StringUtils;

/**
 * Static utility methods that help display results to the user
 *
 * @author Richard Domander
 */
public class ResultUtils {

	private ResultUtils() {}

	/**
	 * Returns a verbal description of the size of the elements in the given
	 * space, e.g. "Area" for 2D images.
	 * 
	 * @param space an N-dimensional space.
	 * @param <S> type of the space.
     * @param <A> type of the axes.
	 * @return the noun for the size of the elements.
	 */
	public static <S extends AnnotatedSpace<A>, A extends TypedAxis> String
		getSizeDescription(S space)
	{
		final long dimensions = AxisUtils.countSpatialDimensions(space);

		if (dimensions == 2) {
			return "Area";
		}
		else if (dimensions == 3) {
			return "Volume";
		}

		return "Size";
	}

	/**
	 * Returns the exponent character of the elements in this space, e.g. '³' for
	 * a spatial 3D space.
	 *
	 * @param space an N-dimensional space.
	 * @param <S> type of the space.
     * @param <A> type of the axes.
	 * @return the exponent character if the space has 2 - 9 spatial dimensions.
	 *         An empty character otherwise.
	 */
	public static <S extends AnnotatedSpace<A>, A extends TypedAxis> char
		getExponent(final S space)
	{
		final long dimensions = AxisUtils.countSpatialDimensions(space);

		if (dimensions == 2) {
			return '\u00B2';
		}
		else if (dimensions == 3) {
			return '\u00B3';
		}
		else if (dimensions == 4) {
			return '\u2074';
		}
		else if (dimensions == 5) {
			return '\u2075';
		}
		else if (dimensions == 6) {
			return '\u2076';
		}
		else if (dimensions == 7) {
			return '\u2077';
		}
		else if (dimensions == 8) {
			return '\u2078';
		}
		else if (dimensions == 9) {
			return '\u2079';
		}

		// Return an "empty" character
		return '\u0000';
	}

	/**
	 * Returns the common unit string that describes the elements in the space.
	 *
	 * @see ResultUtils#getUnitHeader(AnnotatedSpace, UnitService, char)
	 * @param space an N-dimensional space.
	 * @param <S> type of the space.
	 * @param unitService an {@link UnitService} to convert axis calibrations.
	 * @return the unit string with the exponent.
	 */
	public static <S extends AnnotatedSpace<CalibratedAxis>> String getUnitHeader(
		final S space, final UnitService unitService)
	{
		return getUnitHeader(space, unitService, '\u0000');
	}

	/**
	 * Returns the common unit string, e.g. "mm<sup>3</sup>" that describes the
	 * elements in the space.
	 * <p>
	 * The common unit is the unit of the first spatial axis if it can be
	 * converted to the units of the other axes.
	 * </p>
	 *
	 * @param space an N-dimensional space.
	 * @param <S> type of the space.
	 * @param unitService an {@link UnitService} to convert axis calibrations.
	 * @param exponent an exponent to be added to the unit, e.g. '³'.
	 * @return the unit string with the exponent.
	 */
	public static <S extends AnnotatedSpace<CalibratedAxis>> String getUnitHeader(
		final S space, final UnitService unitService, final char exponent)
	{
		final Optional<String> unit = AxisUtils.getSpatialUnit(space, unitService);
		if (!unit.isPresent()) {
			return "";
		}

		final String unitHeader = unit.get();
		if ("pixel".equalsIgnoreCase(unitHeader) || "unit".equalsIgnoreCase(
			unitHeader) || unitHeader.isEmpty())
		{
			// Don't show default units
			return "";
		}

		return "(" + unitHeader + exponent + ")";
	}

	/**
	 * Creates a column for a {@link net.imagej.table.GenericTable} that repeats
	 * the given label on each row.
	 *
	 * @param label the string displayed on each row.
	 * @param rows number of rows created.
	 * @return a column that repeats the label.
	 */
	public static GenericColumn createLabelColumn(final String label, int rows) {
		final GenericColumn labelColumn = new GenericColumn("Label");
		final int n = Math.max(0, rows);
		final String s = StringUtils.isNullOrEmpty(label) ? "-" : label;
		Stream.generate(() -> s).limit(n).forEach(labelColumn::add);
		return labelColumn;
	}

	/**
	 * Creates columns for a {@link net.imagej.table.Table} that describe the
	 * positions of the subspaces in a hyperspace
	 * <p>
	 * For example, if you've split a {X, Y, Z, C, T} space into {X, Y, Z}, the
	 * method returns "Channel" and "Time" columns that list the positions of the
	 * subspaces in C and T.
	 * </p>
	 * 
	 * @see Subspace
	 * @param subspaces the subspaces of a hyperspace.
     * @param <T> type of the elements in the spaces.
	 * @return columns that list the positions of the subspaces.
	 */
	public static <T extends RealType<T> & NativeType<T>> List<LongColumn>
		createCoordinateColumns(List<Subspace<T>> subspaces)
	{
		final List<LongColumn> coordinateColumns = new ArrayList<>();
		if (subspaces == null) {
			return coordinateColumns;
		}
		final AxisType[] types = subspaces.get(0).getAxisTypes().toArray(
			AxisType[]::new);
		final List<long[]> positions = subspaces.stream().map(s -> s.getPosition()
			.toArray()).collect(Collectors.toList());
		for (int i = 0; i < types.length; i++) {
			final AxisType type = types[i];
			final LongColumn coordinateColumn = new LongColumn(type.getLabel());
			final int index = i;
			positions.stream().mapToLong(p -> toConventionalIndex(type, p[index]))
				.forEach(coordinateColumn::add);
			coordinateColumns.add(coordinateColumn);
		}
		return coordinateColumns;
	}

	/**
	 * If needed, converts the given index to the ImageJ1 convention where Z,
	 * Channel and Time axes start from 1.
	 *
	 * @param type type of the axis's dimension.
	 * @param index the index in the axis.
	 * @return index + 1 if type is Z, Channel or Time. Index otherwise.
	 */
	public static long toConventionalIndex(final AxisType type,
		final long index)
	{
		final Stream<AxisType> oneAxes = Stream.of(Axes.Z, Axes.CHANNEL, Axes.TIME);
		if (oneAxes.anyMatch(t -> t.equals(type))) {
			return index + 1;
		}
		return index;
	}
}
