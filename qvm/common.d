module sekitk.qvm.common;

import std.traits: isIntegral, isUnsigned, Signed;
import std.typecons: Tuple;

/*************************************************************
 * Types
 *************************************************************/
alias TypeOfSize= ubyte;
alias TypeOfIndex= ushort;
static assert(8u*ushort.sizeof < (8u*ubyte.sizeof)^^2);

// Enumerate types
enum MatrixType{
  zero= 0u,
  permutation,
  dense,
  band1,
  band3,
  upperTri,
  lowerTri}

/*********************************************
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

/*********************************************
 * Algorithm types of LU decomposition
 *********************************************/
enum DecompScheme: ubyte{
  luWithNonPivDoolittle,			// non-pivoting LU decomposition with Doolittle algorithm
    luWithNonPivCrout,					// non-pivoting LU decomposition with Crout algorithm
    luWithPartialPivDoolittle,	// partial pivoting LU decomposition with Doolittle algorithm
    luWithPartialPivCrout,			// partial pivoting LU decomposition with Crout algorithm
    luWithFullPivDoolittle,			// full pivoting LU decomposition with Doolittle algorithm
    luWithFullPivCrout,					// full pivoting LU decomposition with Crout algorithm
    ldu,				// LDU decomposition
    Cholesky,		// Cholesky decomposition
    modCholesky,	//
    qr,
    singularValue
}

enum DecomposedMat{
  diagonal, upperTri, lowerTri,
  unitary, unitaryLeft, unitaryRight,
  permutationLeft, permutationRight
}

/*********************************************
 * Index type
 *
 * index row i, column j of a_{ij}
 *********************************************/
struct MatrixPosition(TypeOfSize Row, TypeOfSize Column){
  import std.typecons: Tuple;

  invariant(_idxRow <= Row);
  invariant(_idxColumn <= Column);

@safe pure nothrow @nogc:

/***************************
 * constructor
 *
 * Params:
 *	idxRow= 1-based row index
 *	idxColumn= 1-based column index
 ***************************/
this(IdxType)(in IdxType idxRow, in IdxType idxColumn)
if(isIntegral!IdxType)
in{
  assert(idxRow > 0);
  assert(idxColumn > 0);
  assert(idxRow <= Row, "too large row index");
  assert(idxColumn <= Column, "too large column index");
}
do{
  static if(is(typeof(idxRow): TypeOfSize)){
    _idxRow= idxRow;
    _idxColumn= idxColumn;
  }
  else{
    _idxRow= cast(TypeOfSize)(idxRow);
    _idxColumn= cast(TypeOfSize)(idxColumn);
  }
}

  //@disable this();
  this(){
    _idxRow= 1;
    _idxColumn= 1;
  }

  /// copy constructor
  this(ref return scope inout typeof(this) other){}

  const{
    /*************************
     * row index
     *
     * Returns:
     *	1-based row index
     *************************/
    TypeOfSize i() @property{return _idxRow;}

    /*************************
     * column index
     *
     *
     *************************/
    TypeOfSize j() @property{return _idxColumn;}

    /*************************
     *
     *************************/
    Tuple!(TypeOfSize, "i", TypeOfSize, "j") tupleOf(){
      typeof(return) result;
      result.i= _idxRow;
      result.j= _idxColumn;
      return result;
    }

    /*************************
     * explicitly range checking
     *
     * Returns:
     *	true if errors are nothing
     *
     * Throws:
     *	RangeError
     *************************/
    bool rangeCheck(){
      import core.exception: onRangeError;
      typeof(return) result= false;

      assert(_idxRow > 0, "case 1");
      assert(_idxColumn > 0, "case 2");
      assert(_idxRow <= Row, "case 3");
      assert(_idxColumn <= Column, "case 4");
      if(_idxRow > 0 && _idxRow <= Row &&
	 _idxColumn > 0 && _idxColumn <= Column){
	result= true;
      }
      else{
	//throw new InvalidRangeIndex!(row, column)(idxRow, idxColumn);
	onRangeError;
      }
      return result;
    }
  }

package:
  TypeOfSize _idxRow, _idxColumn;
}

/*************************************************************
 * Traits
 *************************************************************/
/*********************************************
 * decomposition
 *********************************************/
template isLU(DecompScheme alg){
  static if(alg is DecompScheme.luWithNonPivDoolittle ||
	    alg is DecompScheme.luWithNonPivCrout ||
	    alg is DecompScheme.luWithPartialPivDoolittle ||
	    alg is DecompScheme.luWithPartialPivCrout ||
	    alg is DecompScheme.luWithFullPivDoolittle ||
	    alg is DecompScheme.luWithFullPivCrout ||
	    alg is DecompScheme.ldu ||
	    alg is DecompScheme.Cholesky ||
	    alg is DecompScheme.modCholesky){
    enum bool isLU= true;
  }
  else{
    enum bool isLU= false;
  }
}

/*********************************************
 * for operator overloading
 *********************************************/
template isPlusOrMinusSign(string Op){
  static if(Op == "+" || Op == "-"){
    enum bool isPlusOrMinusSign= true;
  }
  else{
    enum bool isPlusOrMinusSign= false;
  }
}

/*********************************************
 * Template constraints for each matrix types
 *********************************************/
template matrixConstraint(TypeOfSize Row, TypeOfSize Column, MatrixType Shape){
  static if(Shape is MatrixType.zero || Shape is MatrixType.dense){
    enum bool matrixConstraint= Row > 0 && Column > 0;
  }
  else static if(Shape is MatrixType.permutation ||
		 Shape is MatrixType.band1 ||
		 Shape is MatrixType.upperTri ||
		 Shape is MatrixType.lowerTri){
    enum bool matrixConstraint= Row == Column && Row > 1;
  }
  else static if(Shape is MatrixType.band3){
    enum bool matrixConstraint= Row == Column && Row > 2;
  }
  else{
    static assert(false);
  }
}

/*********************************************
 * Matrix type of its transposed one
 *********************************************/
template ReturnTypeOfTranspose(MatrixType Src){
  static if(Src is MatrixType.upperTri){
    enum MatrixType ReturnTypeOfTranspose= MatrixType.lowerTri;
  }
  else static if(Src is MatrixType.lowerTri){
    enum MatrixType ReturnTypeOfTranspose= MatrixType.upperTri;
  }
  else{
    enum MatrixType ReturnTypeOfTranspose= Src;
  }
}

/*********************************************
 * Return type of matrix multiplications
 *
 * Matrix!(Row, Column, ShapeL) * Matrix!(Column, ColumnR, ShapeR)
 *********************************************/
template ReturnTypeOfOperation(string Op, MatrixType ShapeL, MatrixType ShapeR)
if(Op == "+" || Op == "-"){
  static if(ShapeL == ShapeR){
    enum MatrixType ReturnTypeOfOperation= ShapeL;
  }
  else static if(ShapeL is MatrixType.Zero){
    enum MatrixType ReturnTypeOfOperation= ShapeR;
  }
  else static if(ShapeR is MatrixType.Zero){
    enum MatrixType ReturnTypeOfOperation= ShapeL;
  }
  else{
    static if((ShapeL is MatrixType.band3 && ShapeR is MatrixType.band1) ||
	      (ShapeL is MatrixType.band1 && ShapeR is MatrixType.band3)){
      enum MatrixType ReturnTypeOfOperation= MatrixType.band3;
    }
    else static if((ShapeL is MatrixType.upperTri && ShapeR is MatrixType.band1) ||
		   (ShapeL is MatrixType.band1 && ShapeR is MatrixType.upperTri)){
      enum MatrixType ReturnTypeOfOperation= MatrixType.upperTri;
    }
    else static if((ShapeL is MatrixType.lowerTri && ShapeR is MatrixType.band1) ||
		   (ShapeL is MatrixType.band1 && ShapeR is MatrixType.lowerTri)){
      enum MatrixType ReturnTypeOfOperation= MatrixType.lowerTri;
    }
    else{
      enum MatrixType ReturnTypeOfOperation= MatrixType.dense;
    }
  }
}

/// ditto
template ReturnTypeOfOperation(string OP: "*", MatrixType ShapeL, MatrixType ShapeR){
  static if(ShapeL is MatrixType.zero || ShapeR is MatrixType.zero){
    enum MatrixType ReturnTypeOfOperation= MatrixType.zero;
  }
  else static if(ShapeL is MatrixType.permutation && ShapeR is MatrixType.permutation){
    enum MatrixType ReturnTypeOfOperation= MatrixType.permutation;
  }
  else static if(ShapeL is MatrixType.dense){
    enum MatrixType ReturnTypeOfOperation= ShapeL;
  }
  else static if(ShapeL is MatrixType.band1){
    enum MatrixType ReturnTypeOfOperation= ShapeR;
  }
  else static if(ShapeL is MatrixType.band3){
    static if(ShapeR is MatrixType.band1){
      enum MatrixType ReturnTypeOfOperation= ShapeL;
    }
    else{
      enum MatrixType ReturnTypeOfOperation= MatrixType.dense;
    }
  }
  else static if(ShapeL == MatrixType.upperTri){
    static if(ShapeR is MatrixType.band1 || ShapeR is MatrixType.upperTri){
      enum MatrixType ReturnTypeOfOperation= ShapeL;
    }
    else{
      enum MatrixType ReturnTypeOfOperation= MatrixType.dense;
    }
  }
  else static if(ShapeL == MatrixType.lowerTri){
    static if(ShapeR is MatrixType.band1 || ShapeR is MatrixType.lowerTri){
      enum MatrixType ReturnTypeOfOperation= ShapeL;
    }
    else{
      enum MatrixType ReturnTypeOfOperation= MatrixType.dense;
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
 * Functions
 *************************************************************/
/*********************************************
 *
 * Params:
 *	idxRow= 0-based row index
 *	idxColumn= 0-based column index
 *********************************************/
bool isBijective(TypeOfSize Row,
		 TypeOfSize Column,
		 MatrixType Shape)(in TypeOfSize idxRow, in TypeOfSize idxColumn)
in{
  assert(idxRow >= 0, "too small row index.");
  assert(idxRow < Row, "too large row index.");
  assert(idxColumn >= 0, "too small column index.");
  assert(idxColumn < Column, "too large column index.");
}
do{
  typeof(return) result;

  final switch(Shape){
  case MatrixType.zero, MatrixType.permutation:
    result= false;
    break;
  case MatrixType.dense:
    result= true;
    break;
  case MatrixType.band1:
    result= (idxRow == idxColumn)? true: false;
    break;
  case MatrixType.band3:
    if((idxRow == 0 && idxColumn == 0) ||
       (idxRow == Row-1 && idxColumn == Column-1) ||
       (idxRow+1 == idxColumn) ||
       (idxRow == idxColumn) ||
       (idxRow == idxColumn+1)){
      result= true;
    }
    else{
      result= false;
    }
    break;
  case MatrixType.upperTri:
    if(idxRow <= idxColumn){
      result= true;
    }
    else{
      result= false;
    }
    break;
  case MatrixType.lowerTri:
    if(idxRow >= idxColumn){
      result= true;
    }
    else{
      result= false;
    }
  }
  return result;
}

/*********************************************
 * internal array
 *
 * Params:
 * idxRow= 0-based row index
 * idxColumn= 0-based column index
 *********************************************/
TypeOfIndex indexMap(TypeOfSize Row, TypeOfSize Column,
		     MatrixType Shape,
		     MajorOrder MatOdr)(in uint idxRow,
					in uint idxColumn) @safe pure nothrow @nogc
if(Shape !is MatrixType.zero &&
   Shape !is MatrixType.permutation &&
   matrixConstraint!(Row, Column, Shape))
in{
  assert(idxRow >= 0);
  assert(idxRow < Row, "Row index is invalid.");
  assert(idxColumn >= 0);
  assert(idxColumn < Column, "Column index is invalid.");
  assert(isBijective!(Row, Column, Shape)(cast(TypeOfSize)(idxRow),
					  cast(TypeOfSize)(idxColumn)));
}
do{
  import sekitk.integers.progression: sumFromZero;
  uint result= void;

  static if(Shape is MatrixType.dense){
    final switch(MatOdr){
    case MajorOrder.row:
      result= idxRow*Column+idxColumn;
      break;
    case MajorOrder.column:
      result= idxColumn*Row+idxRow;
    }
  }
  else static if(Shape is MatrixType.band1){
    result= idxRow;
  }
  else static if(Shape is MatrixType.band3){
    result= idxRow*2u+idxColumn;
  }
  else static if(Shape is MatrixType.upperTri){
    final switch(MatOdr){
    case MajorOrder.row:
      result= idxRow*Column+idxColumn-sumFromZero(idxRow);
      break;
    case MajorOrder.column:
      result= sumFromZero(idxRow)+idxColumn;
    }
  }
  else static if(Shape is MatrixType.lowerTri){
    final switch(MatOdr){
    case MajorOrder.row:
      result= sumFromZero(idxRow)+idxColumn;
      break;
    case MajorOrder.column:
      result= idxRow*Column+idxColumn-sumFromZero(idxRow);
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
 * Params:
 *	idxs= MatrixPosition
 *
 * Returns:
 *  [0] index of the specified element is on the internal array
 *  [1] specified position is always zero or not
 *
 * Note:
 *  ROW and COLUMN are assumed less than 0x0100
 ******************************************/
Tuple!(TypeOfIndex, "index",
       bool, "isZero") internalIndexOf(TypeOfSize Row,
				       TypeOfSize Column,
				       MatrixType Shape,
				       MajorOrder MatOdr)(in MatrixPosition!(Row, Column) idxs)
if(Shape !is MatrixType.zero &&
   Shape !is MatrixType.permutation &&
   matrixConstraint!(Row, Column, Shape))
in(idxs.rangeCheck){
  typeof(return) result;

  if(isBijective!(Row, Column, Shape)(cast(TypeOfSize)(idxs.i-1u),
				      cast(TypeOfSize)(idxs.j-1u))){
    result.isZero= false;
    result.index= indexMap!(Row, Column, Shape, MatOdr)(idxs.i-1u, idxs.j-1u);
  }
  else{
    result.isZero= true;
  }

  return result;
}


/*********************************************
 * Kronecker delta
 *********************************************/
Q deltaK(Q)(in Q i, in Q j) @safe pure nothrow @nogc
if(isIntegral!Q){
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
