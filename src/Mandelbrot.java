/*
 * Note that this program uses the Java Vector API, which is only
 * available as an incubator feature in Java 16+. In order to compile
 * and run the program you need to use flag '--add-modules
 * jdk.incubator.vector'. For example:
 *
 * > javac --add-modules jdk.incubator.vector Mandelbrot.java
 *
 * In order to run any program that uses the Mandelbrot class, you
 * similarly have to add the same flag to the java command, e.g.:
 *
 * > java --add-modules jdk.incubator.vector MandelbrotTester
 */

import jdk.incubator.vector.*;
import java.util.Arrays;


// A class for computing the Mandelbrot set and escape times for
// points in the complex plane
public class Mandelbrot {

	// the maximum number of iterations of the system to consider
	private float maxIter = 100.0F;

	// the squared distance to the origin for escape
	private float maxSquareModulus = 100.0F;

	// coordinates of the region to be computed
	private float xMin, xMax, yMin, yMax;

	static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_PREFERRED;

	public Mandelbrot() {

	}

	public Mandelbrot(float maxIter, float maxSquareModulus) {
		this.maxIter = maxIter;
		this.maxSquareModulus = maxSquareModulus;
	}

	public Mandelbrot(float[] params) {
		setAll(params);
	}

	public float getMaxIter() {
		return maxIter;
	}

	// set the region to be considered from points
	public void setRegion(float xMin, float xMax, float yMin, float yMax) {
		this.xMin = xMin;
		this.xMax = xMax;
		this.yMin = yMin;
		this.yMax = yMax;
	}

	// set the region to be considered from an array of coordinates
	public void setRegion(float[] coords) {
		this.xMin = coords[0];
		this.xMax = coords[1];
		this.yMin = coords[2];
		this.yMax = coords[3];
	}

	public void setIterAndModulus(float maxIter, float maxSquaredModulus) {
		this.maxIter = maxIter;
		this.maxSquareModulus = maxSquareModulus;
	}

	// set all parameters; the first four values in params are
	// interpreted as coordinates of the region, while params[4] and
	// params[5] are the maximum iterations and maximum squared
	// modulus, respectively
	public void setAll(float[] params) {
		setRegion(params);
		setIterAndModulus(params[4], params[5]);
	}

	// a baseline implementation of computing escape times for the
	// current region
	// esc is a 2d array to record the escape times,
	// where the first index records rows of the region, and the
	// second index is the column number
	public void escapeTimesBaseline(float[][] esc) {
		float xStep = (xMax - xMin) / esc[0].length;
		float yStep = (yMax - yMin) / esc.length;

		for (int i = 0; i < esc.length; i++) {
			for (int j = 0; j < esc[0].length; j++) {
				int iter = 0;
				float cx = xMin + j * xStep;
				float cy = yMin + i * yStep;

				float zx = 0;
				float zy = 0;

				while (iter < maxIter && zx * zx + zy * zy < maxSquareModulus) {
					float z = zx * zx - zy * zy + cx;
					zy = 2 * zx * zy + cy;
					zx = z;

					iter++;
				}

				esc[i][j] = iter;
			}
		}

	}

	// an optimized implementation of escapeTimesBaseline that uses
	// vector operations
	public void escapeTimesOptimized(float[][] esc) {

		//get distances between each point
		float xStep = (xMax - xMin) / esc[0].length;
		float yStep = (yMax - yMin) / esc.length;


		//initialize an array of steps away from value j
		float[] stepArray = new float[SPECIES.length()];
		stepArray[0] = 0;

		for (int i = 1; i < stepArray.length; i++) {
			stepArray[i] = stepArray[i - 1] + xStep;
		}

		//turn this step array into a vector
		var stepVector = FloatVector.fromArray(SPECIES, stepArray, 0);


		//iterate through each column
		for (int i = 0; i < esc.length; i++) {
			int j = 0;
			//iterate through every set of columns. The number of columns in each batch depends on the number of SPECIES that our computer prefers
			for (; j < SPECIES.loopBound(esc.length); j += SPECIES.length()) {

				//initialize variables into floatVectors.
				var iter = FloatVector.zero(SPECIES);
				var cx = FloatVector.broadcast(SPECIES, xMin + j * xStep).add(stepVector);

				var cy = FloatVector.broadcast(SPECIES, yMin + i * yStep);

				var zx = FloatVector.zero(SPECIES);
				var zy = FloatVector.zero(SPECIES);

				//initialize new floatVectors representing the modulus of a point (x,y)
				var zxsq = zx.mul(zx);
				var zysq = zy.mul(zy);
				var MSModulus = zxsq.add(zysq);

				//create a mask for iterations. This masks checks whether the current iter within a lane is less than maxIter
				var iterMask = iter.compare(VectorOperators.LT, maxIter);

				//create mask for molulus. This mask checks whether the current molulus within a lane is less than maxSquareModulus
				var MSModulusMask = MSModulus.compare(VectorOperators.LT, maxSquareModulus);

				//create mask for the logical intersection between iterMask and MSModulusMask. We will use it to check our conditions
				var executeMask = iterMask.and(MSModulusMask);

				//keep going until none of the lanes satisfy conditions of iter < maxIter and modulus < maxSquareModulus
				while (executeMask.anyTrue()) {

					//calculate modulus
					zxsq = zx.mul(zx);
					zysq = zy.mul(zy);
					MSModulus = zxsq.add(zysq);


					var z = cx.add(zxsq.sub(zysq)); //zx * zx - zy * zy + cx
					zy = cy.add(zx.mul(zy.mul(2.0F))); //2 * zx * zy + cy;
					zx = z;


					//check whether number of iterations still satisfies iter < maxIter and modulus < maxSquareModulus
					iterMask = iter.compare(VectorOperators.LT, maxIter);
					MSModulusMask = MSModulus.compare(VectorOperators.LT, maxSquareModulus);

					//check whether all conditions are filled together
					executeMask = iterMask.and(MSModulusMask);

					//increment iterations
					iter = iter.add(1.0F, executeMask);


				}

				//put iterations into the escape times array
				iter.intoArray(esc[i], j);

			}

		}
	}
}
