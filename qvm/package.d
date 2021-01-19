/*****************************************************************************
 * package "SekiTK"
 *
 * This module provides the linear algebraic data types of quaternion, vector and some matrix.
 * The types are shown in the following list.
 * ---
 * Vector!size                 // vector
 * Vector2                     // 2-dimensional vector
 * Vector3                     // 3-dimensional vector
 * Quaternio                   // quaternion
 * Matrix!(row, column, Shape, MatOdr)  // matrix
 * ---
 *
 * Version: 1.0
 * License: Boost License 1.0
 * Authors: 新井 浩平 (Kohei ARAI), arawi_kohei_takasaki@yahoo.co.jp
 * Source: https://github.com/ARAWI-Kohei-tgj/sekitk/src/qvm/package.d
 *****************************************************************************/
module sekitk.qvm;

import std.traits: isFloatingPoint;
import sekitk.complex.pseudo: isComplex, BaseRealType;
public import sekitk.qvm.common: MatrixType, MajorOrder, DecompScheme, DecomposedMat, MatrixPosition, arrayLength;
public import sekitk.qvm.exception: ZeroNorm, SingularMatrix;
//public import sekitk.algorithm.decomposition;	// can be used in standalone

/*************************************************************
 * Template SekiTK!(T, T threshold)
 *
 *
 *	T = floating point type in $(B R) or $(B C)
 *	Threshold = a number which its magnitude is less than "threshold" is treated as zero
 *************************************************************/
template SekiTK(T, real Threshold= 1.0e-6L)
if(isFloatingPoint!T || isComplex!T){
	// Data types
	import sekitk.qvm.euclideanvector;
	import sekitk.qvm.matrix;
	import sekitk.qvm.quaternion;

	mixin VectorImpl!(T, Threshold);	// Vector
	mixin MatrixImpl!(T, Threshold);	// Matrix
	static if(isFloatingPoint!T){
		import sekitk.qvm.quaternion;
		alias Vector2= Vector!2u;	// Euclidean vector R^2
		alias Vector3= Vector!3u;	// Euclidean vector R^3
		mixin QuaternionImpl!(T, Threshold);	/// Quaternion
	}

	// Functions
	@safe pure nothrow{
		/******************************************
		 * Identity matrix (any a[i=j]= 1, others= 0)
		 *
		 * In order to decrease the memory consumption,
		 * a n-sized identity matrix should be constructed
		 * as diagonal matrix.
		 * However, if you want to use it with other
		 * type of matrix, this function will be useful.
		 ******************************************/
		Matrx!(Size, Size, Shape, MatOdr) identityMatrix(TypeOfSize Size,
																										 MatrixType Shape,
																										 MajorOrder MatOdr)(){
			import sekitk.qvm.common: arrayLength;
		  import sekitk.qvm.indexset: IndexSetDiag, IndexSetStrictTriR, IndexSetStrictTriL;
			T[arrayLength!(Size, Size, Shape)] buf= void;
			static if(Shape is MatrixType.diag){
				buf[]= T(1.0L);
			}
			else{
				{
					auto idxSet0= IndexSetStrictTriR!(Size, Shape, MatOdr)();
					auto idxSet1= IndexSetStrictTriL!(Size, Shape, MatOdr)();
					foreach(TypeOfSize idx; merge(idxSet0, idxSet1)) buf[idx]= T(0.0L);
				}
				{
					auto idxSet= IndexSetDiag!(Size, Size, Shape, MatOdr)();
					foreach(TypeOfSize idx; idxSet) buf[idx]= T(1.0L);
				}
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
			const mat= new Matrix!(3, 3, MatrixType.dense, MajorOrder.column)(temp);
			return mat.det;
		}
	}

	/**************************
	 * vector triple product for C^3
	 *
	 * a \tiems (b \times c)
	 **************************/
	Vector!3 vectorTripleProd(in Vector!3 a,
														in Vector!3 b,
														in Vector!3 c) @safe pure nothrow @nogc{
		const T[2] coeff= [a*c, a*b];
		return coeff[0]*b-coeff[1]*c;
	}
}
// end of template `SekiTK'


