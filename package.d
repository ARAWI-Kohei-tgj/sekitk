/*****************************************************************************
 * package "SekiTK"
 *
 * This module contains the definitions of vector, matrix and quaternion.
 * The provided types are shown in the following list.
 * ---
 * Vector!size                 // vector
 * Vector2                     // 2-dimensional vector
 * Vector3                     // 3-dimensional vector
 * Quaternio                   // quaternion
 * Matrix!(row, column)        // rectangural matrix
 * Matrix!size                 // square matrix
 * UpperTriangularMatrix!size  // upper triangural matrix
 * LowerTriangularMatrix!size  // lower triangural matrix
 * BandMatrix!size             // band matrix
 * DiagonalMatrix!size         // diagonal matrix
 * ---
 *
 * Version: 0.0
 *****************************************************************************/
module sekitk;

import std.traits: isFloatingPoint;
import numeric.pseudocmplx: isComplex, BaseRealType;
public import sekitk.base: MatrixType, MajorOrder, DecompScheme, DecomposedMat, MatrixPosition, arrayLength;
public import sekitk.exceptions;
public import sekitk.lu;	// can be used in standalone

/*************************************************************
 * Template SekiTK!(T, T threshold)
 *
 * Params:
 *	T = floating point (R or C)
 *	threshold = a number which its magnitude is less than "threshold" is treated as zero
 *  matodr= row major order or column major order
 *************************************************************/
template SekiTK(T, real Threshold= 1.0e-6)
if(isFloatingPoint!T || isComplex!T){
	/********************************************
	 * Class mixin
	 ********************************************/
	import sekitk.euclideanvector;
	import sekitk.matrix;
	import sekitk.quaternion;

	mixin VectorImpl!(T, Threshold);	/// Vector
	mixin MatrixImpl!(T, Threshold);	/// Matrix
	static if(isFloatingPoint!T){
		import sekitk.quaternion;
		alias Vector2= Vector!2u;	// Euclidean vector R^2
		alias Vector3= Vector!3u;	// Euclidean vector R^3
		mixin QuaternionImpl!(T, Threshold);	/// Quaternion
	}

	/********************************************
	 * Functions
	 ********************************************/
	@safe pure nothrow{
		/**************************
		 * Identity matrix (any a[i=j]= 1, others= 0)
		 **************************/
		Matrx!(Size, Size, Shape) identityMatrix(uint Size, MatrixType Shape)(){
			import numeric.sekitk.base: arrayLength, idxSetDiagImpl;
			T[arrayLength!(Size, Size, Shape)] buf= void;
			static if(Shape is MatrixType.diag){
				buf[]= T(1.0);
			}
			else{
				buf[]= T(0.0);
				auto idxset= new IndexSetDiag!(Size, Size, Shape);
				foreach(idx; idxset) buf[idx]= T(1.0);
			}
			return new typeof(return)(buf);
		}

		/**************************
		 * scalar triple product for C^3
		 *
		 * a \cdot (b \times c)
		 **************************/
		T scalarTripleProd(in Vector!3 a, in Vector!3 b, in Vector!3 c){
			const Vector!3u[3] temp= [a, b, c];
			const auto mat= new SekiTK!(T, Threshold).Matrix!(3, 3, MatrixType.dense, MajorOrder.column)(temp);
			return mat.det;
		}

		/**************************
		 * vector triple product for C^3
		 *
		 * a \tiems (b \times c)
		 **************************/
		Vector!3 vectorTripleProd(in Vector!3 a, in Vector!3 b, in Vector!3 c){
			const T[2] coeff= [a*c, a*b];
			return coeff[0]*b-coeff[1]*c;
		}
	}
}
// end of template `SekiTK'


