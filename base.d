/*****************************************************************************
 * Base module of SekiTK
 *
 * Author: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 *****************************************************************************/
module sekitk.base;

import std.traits: isIntegral, isUnsigned, Signed;


import std.typecons: Tuple;

/**************************************************************
 * Types
 **************************************************************/
alias TypeOfSize= ubyte;
alias TypeOfIndex= ushort;
static assert(8u*ushort.sizeof < (8u*ubyte.sizeof)^^2);

// Enumerate types
/**************************************************************
 * Type of matrix
 *
 * Zero= zero matrix
 * Dense= dense matrix
 * Band1= diagonal matrix
 * Band3= tridiagonal matrix
 * UpperTri= upper triangular matrix
 * LowerTri= lower triangular matrix
 **************************************************************/
enum MatrixType{
	zero= 0u,
	dense,
	band1,
	band3,
	upperTri,
	lowerTri}

/**********************************************
 * Straging type of elements
 *
 * Row= row major order
 *  0  1  2  3
 *  4  5  6  7
 *  8  9 10 11
 * 12 13 14 15
 *
 * Column= column major order
 *  0  4  8 12
 *  1  5  9 13
 *  2  6 10 14
 *  3  7 11 15
 *********************************************/
enum MajorOrder{row, column}

/**********************************************
 * List of matrix decomposition methods
 *
 * LUP -> PA= LU
 *		partial pivoting LU decomposition
 *
 * LUPQ -> PAQ= LU
 *		full pivoting LU decomposition
 *
 * LDU -> A= LDU
 * Choresky
 * QR
 * SV
 *********************************************/
enum DecompScheme{
	singularValue,
	qr,
	lup,
	lupq,
	ldu,
	cholesky}

enum DecomposedMat{
  diagonal, upperTri, lowerTri,
	unitary, unitaryLeft, unitaryRight,
	permutationLeft, permutationRight}

// Constraints
/**********************************************
 * Template constraints for each matrix types
 *********************************************/
template matrixConstraint(TypeOfSize Row, TypeOfSize Column, MatrixType Shape){
	static if(Shape is MatrixType.zero || Shape is MatrixType.dense){
		enum bool matrixConstraint= (Row > 0 && Column > 0)? true : false;
	}
	else static if(Shape is MatrixType.band1 ||
								 Shape is MatrixType.upperTri ||
								 Shape is MatrixType.lowerTri){
		enum bool matrixConstraint= (Row == Column && Row > 1)? true : false;
	}
	else static if(Shape is MatrixType.band3){
		enum bool matrixConstraint= (Row == Column && Row > 2)? true : false;
	}
	else{
		static assert(false);
	}
}

/**************************************************************
 * Matrix type of its transposed one
 **************************************************************/
template MatrixTypeOfTranspose(MatrixType Src){
	static if(Src is MatrixType.upperTri){
		enum MatrixType MatrixTypeOfTranspose= MatrixType.lowerTri;
	}
	else static if(Src is MatrixType.lowerTri){
		enum MatrixType MatrixTypeOfTranspose= MatrixType.upperTri;
	}
	else{
		enum MatrixType MatrixTypeOfTranspose= Src;
	}
}

/*********************************************
 * Return type of matrix multiplications
 *
 * Matrix!(Row, Column, ShapeL) * Matrix!(Column, ColumnR, ShapeR)
 *********************************************/
template MatOpReturnType(string Op, MatrixType ShapeL, MatrixType ShapeR)
if(Op == "+" || Op == "-"){
	static if(ShapeL == ShapeR){
		enum MatrixType MatOpReturnType= ShapeL;
	}
	else static if(ShapeL is MatrixType.Zero){
		enum MatrixType MatOpReturnType= ShapeR;
	}
	else static if(ShapeR is MatrixType.Zero){
		enum MatrixType MatOpReturnType= ShapeL;
	}
	else{
		static if((ShapeL is MatrixType.band3 && ShapeR is MatrixType.band1) ||
							(ShapeL is MatrixType.band1 && ShapeR is MatrixType.band3)){
			enum MatrixType MatOpReturnType= MatrixType.band3;
		}
		else static if((ShapeL is MatrixType.upperTri && ShapeR is MatrixType.band1) ||
									 (ShapeL is MatrixType.band1 && ShapeR is MatrixType.upperTri)){
			enum MatrixType MatOpReturnType= MatrixType.upperTri;
		}
		else static if((ShapeL is MatrixType.lowerTri && ShapeR is MatrixType.band1) ||
									 (ShapeL is MatrixType.band1 && ShapeR is MatrixType.lowerTri)){
			enum MatrixType MatOpReturnType= MatrixType.lowerTri;
		}
		else{
			enum MatrixType MatOpReturnType= MatrixType.dense;
		}
	}
}

/// ditto
template MatOpReturnType(string OP: "*", MatrixType ShapeL, MatrixType ShapeR){
	static if(ShapeL is MatrixType.zero || ShapeR is MatrixType.zero){
		enum MatrixType MatOpReturnType= MatrixType.zero;
	}
	else static if(ShapeL is MatrixType.dense){
		enum MatrixType MatOpReturnType= ShapeL;
	}
	else static if(ShapeL is MatrixType.band1){
		enum MatrixType MatOpReturnType= ShapeR;
	}
	else static if(ShapeL is MatrixType.band3){
		static if(ShapeR is MatrixType.band1){
			enum MatrixType MatOpReturnType= ShapeL;
		}
		else{
			enum MatrixType MatOpReturnType= MatrixType.dense;
		}
	}
	else static if(ShapeL == MatrixType.upperTri){
		static if(ShapeR is MatrixType.band1 || ShapeR is MatrixType.upperTri){
			enum MatrixType MatOpReturnType= ShapeL;
		}
		else{
			enum MatrixType MatOpReturnType= MatrixType.dense;
		}
	}
	else static if(ShapeL == MatrixType.lowerTri){
		static if(ShapeR is MatrixType.band1 || ShapeR is MatrixType.lowerTri){
			enum MatrixType MatOpReturnType= ShapeL;
		}
		else{
			enum MatrixType MatOpReturnType= MatrixType.dense;
		}
	}
}

/*********************************************
 * Length of internal static array
 *
 * Params:
 * 	row = number of rows
 * 	column = number of columns
 * 	shape = type of matrix
 *********************************************/
template arrayLength(TypeOfSize Row, TypeOfSize Column, MatrixType Shape)
if(matrixConstraint!(Row, Column, Shape)){
	static if(Shape is MatrixType.zero){
		enum TypeOfIndex arrayLength= 0u;
	}
	else static if(Shape is MatrixType.dense){
		enum TypeOfIndex arrayLength= Row*Column;
	}
	else static if(Shape is MatrixType.band1){
		enum TypeOfIndex arrayLength= Row;
	}
	else static if(Shape is MatrixType.band3){
		enum TypeOfIndex arrayLength= 3u*(Row-2u)+4u;
	}
	else static if(//Shape is MatrixType.Hermitian ||
								 Shape is MatrixType.upperTri ||
								 Shape is MatrixType.lowerTri
								 ){
		enum TypeOfIndex arrayLength= Row*(Row+1u)/2u;
	}
/+
	else static if(Shape is MatrixType.SkewHermitian){
		static if(isComplex!T){
			enum TypeOfIndex arrayLength= Row*(Row-1u)/2u;
		}
		else{
			enum TypeOfIndex arrayLength= Row*(Row+1u)/2u;
		}
}+/
	else{
		static assert(false);
	}
}
@safe pure nothrow @nogc unittest{
	assert(arrayLength!(4, 5, MatrixType.dense) == 20u);	// rectangular dense matrix
	assert(arrayLength!(3, 3, MatrixType.band1) == 3u);	// diagonal matrix
	assert(arrayLength!(4, 4, MatrixType.band3) == 10u);	// tridiagonal matrix
	assert(arrayLength!(4, 4, MatrixType.upperTri) == 10u);	// upper triangular matrix
	assert(arrayLength!(5, 5, MatrixType.lowerTri) == 15u);	// lower triangular matrix
}

/*************************************************************
 * Index type
 *
 * index row i, column j of a_{ij}
 *************************************************************/
struct MatrixPosition(TypeOfSize row, TypeOfSize column){
	import std.typecons: Tuple;
	/****************************
	 * Constructor
	 ****************************/
	this(in size_t i, in size_t j) @safe pure nothrow @nogc{
		import core.exception: onRangeError;
		if(i > typeof(row).max || j > typeof(column).max){
			onRangeError;
		}
		else{
			idxRow= cast(TypeOfSize)i;
			idxColumn= cast(TypeOfSize)j;
		}
	}

	@safe pure const @property{
		auto i() nothrow @nogc{return idxRow;}
	  auto j() nothrow @nogc{return idxColumn;}

		Tuple!(TypeOfSize, "i", TypeOfSize, "j") tupleOf() nothrow @nogc{
			typeof(return) result;
			result.i= idxRow;
			result.j= idxColumn;
			return result;
		}

		bool rangeCheck() nothrow @nogc{
			//import numeric.sekitk.exceptions: InvalidRangeIndex;
			import core.exception: onRangeError;
			typeof(return) result;
			assert(result == false, "specification change or compiler bug");

			if(idxRow > 0u && idxRow <= row &&
				 idxColumn > 0u && idxColumn <= column){
				result= true;
			}
			else{
				//throw new InvalidRangeIndex!(row, column)(idxRow, idxColumn);
				onRangeError;
			}
			return result;
		}
	}

protected:
	TypeOfSize idxRow, idxColumn;
}

/**
 *
 ****/
bool isBijective(TypeOfSize Row, TypeOfSize Column,
								 MatrixType Shape)(in TypeOfSize idxRow, in TypeOfSize idxColumn)
in(idxRow < Row, "Row index is out of bound.")
in(idxColumn < Column, "Column index is out of bound."){
	typeof(return) result;

	final switch(Shape){
	case MatrixType.zero:
		result= false;
		break;
  case MatrixType.dense:
		result= true;
		break;
  case MatrixType.band1:
		result= (idxRow == idxColumn)? true: false;
		break;
  case MatrixType.band3:
		if((idxRow == 0 && idxColumn == 1) ||
			 (idxRow == Row-1 && idxColumn == Column-2) ||
			 (idxRow == idxColumn-1) ||
			 (idxRow == idxColumn) ||
			 (idxRow == idxColumn+1)){
			result= true;
		}
		else{
			result= false;
		}
		break;
	case MatrixType.upperTri:
		break;
	case MatrixType.lowerTri:
	}
	return result;
}

/****************
 * internal array
 *
 * Params:
 * idxRow= row index (0 ≤ idxRow < Row)
 * idxColumn= column index (0 ≤ idxColumn < Column)
 ***/
TypeOfIndex indexMap(TypeOfSize Row, TypeOfSize Column,
										 MatrixType Shape,
										 MajorOrder MatOdr)(in uint idxRow,
																				in uint idxColumn) @safe pure nothrow @nogc
if(Shape !is MatrixType.zero &&
	 matrixConstraint!(Row, Column, Shape))
in(idxRow < Row, "Row index is invalid.")
in(idxColumn < Column, "Column index is invalid.")
in(isBijective!(Row, Column, Shape)(cast(TypeOfSize)(idxRow),
																		cast(TypeOfSize)(idxColumn))){
	import mathematic.progression: sumFromZero;
	uint result;

	static if(Shape is MatrixType.dense){
		final switch(MatOdr){
		case MajorOrder.row:
			result= cast(TypeOfIndex)(idxRow*Column+idxColumn);
			break;
		case MajorOrder.column:
			result= cast(TypeOfIndex)(idxColumn*Row+idxRow);
			break;
		case MajorOrder.diag:
			uint tempIdx= idxColumn-1;
			const uint en= 1u+ ((idxRow-idxColumn >= 0)?
													idxRow-idxColumn:
													idxColumn-idxRow);
			if(idxRow-idxColumn < 0) tempIdx -= Row;
			foreach(k; 0..en) tempIdx += 2*(Row-k)-1;
			result= cast(TypeOfIndex)tempIdx;
		}
	}
	else static if(Shape is MatrixType.band1){
		result= idxRow;
	}
	else static if(Shape is MatrixType.band3){
		final switch(MatOdr){
		case MajorOrder.row, MajorOrder.column:
			result= idxRow*2u+idxColumn;
			break;
		case MajorOrder.diag:
			if(idxRow == idxColumn) result= idxRow;
			else if(idxRow < idxColumn) result= Row-1 +idxColumn;
			else result= 2*Row-2+idxRow;
		}
	}
	else static if(Shape is MatrixType.upperTri){
		final switch(MatOdr){
		case MajorOrder.row:
			result= cast(TypeOfIndex)(idxRow*Column+idxColumn-sumFromZero(idxRow));
			break;
		case MajorOrder.column:
			result= cast(TypeOfIndex)(sumFromZero(idxRow)+idxColumn);
			break;
		case MajorOrder.diag:
		}
	}
	else static if(Shape is MatrixType.lowerTri){
		final switch(MatOdr){
		case MajorOrder.row:
			result= cast(TypeOfIndex)(sumFromZero(idxRow)+idxColumn);
			break;
		case MajorOrder.column:
			result= cast(TypeOfIndex)(idxRow*Column+idxColumn-sumFromZero(idxRow));
			break;
		case MajorOrder.diag:
		}
	}
	else{
		assert(false);
	}

	assert(result <= TypeOfIndex.max);
	return cast(TypeOfIndex)result;
}

/******************************************
 * Implementation for random access to internal array 
 * specified by index[i, j]
 * 
 * Returns:
 *  [0] index of the specified element on the internal array
 *  [1] specified position is always zero or not
 *
 * Note:
 *  ROW and COLUMN are assumed less than 0x0100
 ******************************************/
 Tuple!(TypeOfIndex, "index", bool, "isZero") indexOfInternalArray(TypeOfSize Row, TypeOfSize Column, MatrixType Shape, MajorOrder MatOdr)(in MatrixPosition!(Row, Column) idxs)
if(matrixConstraint!(Row, Column, Shape))
in(idxs.rangeCheck){
	const TypeOfSize i= cast(ubyte)(idxs.i-1), j= cast(ubyte)(idxs.j-1);
	typeof(return) result;

	if(isBijective!(Row, Column, Shape)(i, j)){
		result.isZero= false;
		result.index= indexMap!(Row, Column, Shape, MatOdr)(i, j);
	}
	else{
		result.isZero= true;
	}

	return result;
}

/**
 *
 **/


/**
 * .idup
 */
auto trustedAssumeUnique(U)(U t) @trusted pure nothrow @nogc{
	import std.exception: assumeUnique;
	return assumeUnique(t);
}

/******************************
 * for operator overloading
 *****************************/
template isPlusOrMinusSign(string Op){
	static if(Op == "+" || Op == "-") enum bool isPlusOrMinusSign= true;
	else enum bool isPlusOrMinusSign= false;
}

/**
 *
 **/
template isUnsignedInt(T){
	import std.traits: isIntegral, isUnsigned;

	static if(isIntegral!T && isUnsigned!T){
		enum bool isUnsignedInt= true;
	}
	else{
		enum bool isUnsignedInt= false;
	}
}
@safe pure nothrow @nogc unittest{
	assert(isUnsignedInt!ushort);
	assert(!isUnsignedInt!long);
	assert(!isUnsignedInt!double);
}

template isSignedInt(T){
	import std.traits: isIntegral, isUnsigned;

	static if(isIntegral!T && !isUnsigned!T){
		enum bool isSignedInt= true;
	}
	else{
		enum bool isSignedInt= false;
	}
}
@safe pure nothrow @nogc unittest{
	assert(!isSignedInt!ushort);
	assert(isSignedInt!long);
	assert(!isSignedInt!double);
}

/**
 * Kronecker delta
 */
T deltaK(T)(in T i, in T j) @safe pure nothrow @nogc
if(isIntegral!T && isUnsigned!T){
	return (i == j)? 1u : 0u;
}

/**
 * Levi-Civita symbol
 **/
/+
Signed!T epsilonE(T, uint ODR)(in U i, in U j) @safe pure nothrow @nogc
if(isIntegral!T && isUnsigned!T && ODR >= 2u){
	static if(ODR == 2u){}
}+/
