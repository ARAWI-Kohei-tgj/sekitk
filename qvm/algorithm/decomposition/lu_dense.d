module sekitk.qvm.algorithm.decomposition.lu_dense;

/**************************************************************
 * LU decomposition
 *
 * This template is a dummy for N=1 dense matrix.
 **************************************************************/
enum string LU_DENSE_IMPL1= q{
  /********************************************
   * A dummy function for the decomposition
   ********************************************/
  Result!Method decompose(DecompScheme Method,
			  ArgType)(in ArgType mat) @safe pure nothrow @nogc
  if(is(ArgType: T) || is(ArgType: T[1])){
    import std.traits: isArray;
    static if(isArray!ArgType){
      return new typeof(return)(mat[0]);
    }
    else{
      return new typeof(return)(mat);
    }
  }

  /********************************************
   * Data type for result of the decomposition
   ********************************************/
  struct Result(DecompScheme Method){
    @safe pure:
    // Constructors
    nothrow @nogc{
      this(in T mat){
        lu= mat;
      }

      this(ref return scope inout typeof(this) other){}
    }

    // Methods
    const{
      /************************
       * Returns:
       *	true if invertible
       ************************/
      @property bool isInvertible() nothrow @nogc{
        return approxEqualZero!(T, Threshold)(lu)? false: true;
      }

      /// determinant
      alias det= lu;

      /************************
       * Inverse
       *
       * Throws:
       *	SingularMatrix
       ************************/
      T inverseArray() @property{
	typeof(return) result= Identity!(T, "*")/lu;
	if(isNaN(result)){
	  throw new SingularMatrix();
	}
	return result;
      }
    }

  private:
    T lu;
  }
};

/**************************************************************
 * LU decomposition
 *
 * for dense matrix (N > 1)
 **************************************************************/
enum string LU_DENSE_IMPL2= q{
  import std.algorithm: merge;
  import std.array: staticArray;
  import std.range: iota, zip;

  /********************************************
   * This function factrizes the square matrix A to L, U and P using Doolittle algorithm
   *
   * $(I $(B A))= $(I $(B LU))
   * $(I $(B PA))= $(I $(B LU))
   * $(I $(B PAQ))= $(I $(B LU))
   ********************************************/
  Result!Method decompose(DecompScheme Method,
			  ArgType)(in ArgType mat) @safe pure nothrow @nogc
  if(isNormalLU!Method && (is(ArgType: T[arrayLength!(Size, Size, Shape)]) ||
			   is(ArgType: T[Size][Size]))){
    import std.algorithm: maxIndex, swap, swapRanges;
    import std.math: abs;

    static if(is(ArgType: T[arrayLength!(Size, Size, Shape)])){
      T[Size][Size] a= arrayTo2Dim!(T, Size, Size)(mat);
    }
    else{
      T[Size][Size] a= mat;
    }

    bool invertible= true;

    static if(isPartialPivotting!Method){
      enum string PIV_IMPL= PARTIAL_PIV_IMPL;

      TypeOfSize[Size] p= iota(Size).staticArray!Size;
      size_t idxMax= void;
      TypeOfSize counterEx= 0;
    }
    else static if(isFullPivotting!Method){
      enum string PIV_IMPL= q{};

      TypeOfSize[Size][2] p= iota(Size).staticArray!Size;
      size_t[2] idxMax= void;
      TypeOfSize[2] counterEx= 0;
    }
    else{
      enum string PIV_IMPL= "";
    }

    foreach(scope TypeOfSize idxD; 0u..Size){	// row major order -> seeks horizontally
      mixin(PIV_IMPL);	// pivotting

      static if(isDoolittle!Method){
	foreach(scope i; idxD+1u..Size){	// row major order -> seek vertically
	  a[i][idxD] /= a[idxD][idxD];
	  foreach(scope k; idxD+1u..Size) a[i][k] -= a[i][idxD]*a[idxD][k];
	}
      }
      else static if(isCrout!Method){
	foreach(scope i; idxD+1u..Size){
	  a[idxD][i] /= a[idxD][idxD];
	  foreach(scope k; idxD+1u..Size) a[k][i] -= a[k][idxD]*a[idxD][i];
	}
      }
      else{
	static assert(false);
      }
    }

    static if(isPartialPivotting!Method || isFullPivotting!Method){
      return typeof(return)(a, invertible, p, counterEx);
    }
    else{
      return typeof(return)(a, invertible);
    }
  }

  /******************************************
   * Data type for result of the decomposition
   ******************************************/
  struct Result(DecompScheme Method)
  if(isPartialPivotting!Method){
    import std.typecons: Tuple, tuple, Ternary;
    @safe pure:
    // constructors
    nothrow @nogc{
      /// copy constructor
      this(ref return scope inout typeof(this) other){}

      /************************
       * Params:
       *	isInvertible= 
       *	mat= 
       *	perm= number of the permutation
       *	ex= exchange number of row or column
       ************************/
      this(in T[Size][Size] mat, in bool isInvertible,
	   in ubyte[Size] perm, in TypeOfSize ex){
	_invertible= isInvertible;
	foreach(scope TypeOfSize i; 0u..Size)
	  foreach(scope TypeOfSize j; 0u..Size){{
	      _values[indexMap!(Size, Size, Shape, MatOdr)(i, j)]= mat[i][j];
	    }}
	_sequence[]= perm[];
	_counterEx= ex;
      }
    }

    // Other methods
    const{
      /**********************
       * 
       **********************/
      bool isInvertible() nothrow @nogc @property{
	return _invertible;
      }

      /************************
       * Determinant
       ************************/
      T det() nothrow @nogc{
	import std.range: dropOne;
	typeof(return) result= _values[0];
	auto idxSet= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
	foreach(scope TypeOfIndex idx; idxSet.dropOne) result *= _values[idx];
	if(_counterEx%2) result= -result;
	return result;
      }

      /************************
       * Inverse
       *
       * Throws:
       *	SingularMatrix
       ************************/
      T[arrayLength!(Size, Size, MatrixType.dense)] inverseArray(){
	import sekitk.qvm.exception: SingularMatrix;
	typeof(return) result= void;

	if(_invertible){
	  T[Size][Size] inv2d;

	  foreach(scope TypeOfSize j; 0u..Size){
	    foreach(scope TypeOfSize i; 0u..Size){
	      inv2d[i][j]= (_sequence[i] == j)? VALUE_ONE: VALUE_ZERO;
	      foreach(scope TypeOfSize k; 0u..i) inv2d[i][j] -= _values[indexMap!(Size, Size, Shape, MatOdr)(i, k)]*inv2d[k][j];
	    }

	    foreach_reverse(scope i; 0u..Size){	// NOTE: it must be an down counting
	      foreach(scope k; i+1u..Size) inv2d[i][j] -= _values[indexMap!(Size, Size, Shape, MatOdr)(i, k)]*inv2d[k][j];
	      inv2d[i][j] /= _values[indexMap!(Size, Size, Shape, MatOdr)(i, i)];
	    }
	  }

	  size_t st;
	  foreach(scope i; 0..Size){
	    st= i*Size;
	    result[st..st+Size]= inv2d[i];
	  }
	}
	else{
	  throw new SingularMatrix();
	}

	return result;
      }

      /**********************
       * Lower triangular matrix $(I L)
       **********************/
      T[arrayLength!(Size, Size, MatrixType.lowerTri)]
	matrix(DecomposedMat Mat: DecomposedMat.lowerTri)() nothrow @nogc{
	enum TypeOfIndex LEN= arrayLength!(Size, Size, MatrixType.lowerTri);
	T[LEN] result= void;

	static if(isDoolittle!Method){
	  {	// diagonal elements
	    auto idxSetDest= IndexSetDiag!(Size, Size, MatrixType.lowerTri, MatOdr)();
	    foreach(scope TypeOfIndex idx; idxSetDest) result[idx]= VALUE_ONE;
	  }
	  {	// strict lower triangular elements
	    auto idxSetSrc= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
	    auto idxSetDest= IndexSetStrictTriL!(Size, MatrixType.lowerTri, MatOdr)();
	    foreach(scope idxSrc, idxDest;
		    zip(idxSetSrc, idxSetDest)) result[idxDest]= _values[idxSrc];
	  }
	}
	else{
	  auto idxSetSrc0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
	  auto idxSetSrc1= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
	  foreach(sope idxSrc, idxDest;
		  zip(merge(idxSetSrc0, idxSetSrc1), iota(LEN))
		  ) result[idxDest]= _values[idxSrc];
	}

	return result;
      }

      /************************
       * Upper triangular matrix $(I U)
       ************************/
      T[arrayLength!(Size, Size, MatrixType.upperTri)]
	matrix(DecomposedMat Mat: DecomposedMat.upperTri)() nothrow @nogc{
	enum TypeOfIndex LEN= arrayLength!(Size, Size, MatrixType.upperTri);
	T[LEN] result= void;

	static if(isDoolittle!Method){
	  auto idxSetSrc0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
	  auto idxSetSrc1= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
	  foreach(scope idxSrc, idxDest;
		  zip(merge(idxSetSrc0, idxSetSrc1), iota(LEN))
		  ) result[idxDest]= _values[idxSrc];
	}
	else{
	  {	// diagonal elements
	    auto idxSetDest= IndexSetDiag!(Size, Size, MatrixType.upperTri, MatOdr)();
	    foreach(scope TypeOfIndex idx; idxSetDest) result[idx]= VALUE_ONE;
	  }
	  {	// strict upper triangular elements
	    auto idxSetSrc= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
	    auto idxSetDest= IndexSetStrictTriR!(Size, MatrixType.upperTri, MatOdr)();
	    foreach(scope idxSrc, idxDest; zip(idxSetSrc, idxSetDest)) result[idxDest]= _values[idxSrc];
	  }
	}

	return result;
      }

      /************************
       * Permutation Matrix $(I P)
       ************************/
      auto matrix(DecomposedMat Mat: DecomposedMat.permutationLeft)() nothrow{
	return tuple!("sequence", "isOdd")(_sequence, Ternary(_counterEx%2 > 0));
      }
    }

  private:
    TypeOfSize[Size] _sequence;
    T[arrayLength!(Size, Size, MatrixType.dense)] _values;
    TypeOfSize _counterEx= 0;
    bool _invertible= false;
  }

  /// ditto
  struct Result(DecompScheme Method: DecompScheme.luWithNonPivDoolittle){
    import std.typecons: Tuple, tuple, Ternary;
    @safe pure:
    // constructors
    nothrow @nogc{
      /// copy constructor
      this(ref return scope inout typeof(this) other){}

      /************************
       * Params:
       *	isInvertible= 
       *	mat= 
       *	perm= number of the permutation
       *	ex= exchange number of row or column
       ************************/
      this(in bool isInvertible, in T[Size][Size] mat){
	_invertible= isInvertible;
	foreach(scope TypeOfSize i; 0u..Size)
	  foreach(scope TypeOfSize j; 0u..Size){{
	      _values[indexMap!(Size, Size, Shape, MatOdr)(i, j)]= mat[i][j];
	    }}
      }
    }

    // Other methods
    const{
      /************************
       * Returns true if the matrix is invertible, otherwise false
       ************************/
      bool isInvertible() nothrow @nogc @property{
	return _invertible;
      }

      /************************
       * Determinant
       ************************/
      T det() nothrow @nogc{
	import std.range: dropOne;
	typeof(return) result= _values[0];
	auto idxSet= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
	foreach(scope TypeOfIndex idx; idxSet.dropOne) result *= _values[idx];
	return result;
      }

      /************************
       * Inverse
       *
       * Throws:
       *	SingularMatrix
       ************************/
      T[arrayLength!(Size, Size, MatrixType.dense)] inverseArray(){
	import sekitk.qvm.exception: SingularMatrix;
	typeof(return) result= void;

	if(_invertible){
	  T[Size][Size] inv2d;

	  foreach(scope TypeOfSize j; 0u..Size){
	    foreach(scope TypeOfSize i; 0u..Size){
	      inv2d[i][j]= (i == j)? VALUE_ONE: VALUE_ZERO;	// FIXME:
	      foreach(scope TypeOfSize k; 0u..i) inv2d[i][j] -= _values[indexMap!(Size, Size, Shape, MatOdr)(i, k)]*inv2d[k][j];
	    }

	    foreach_reverse(scope i; 0u..Size){	// it must be an down counting
	      foreach(scope k; i+1u..Size) inv2d[i][j] -= _values[indexMap!(Size, Size, Shape, MatOdr)(i, k)]*inv2d[k][j];
	      inv2d[i][j] /= _values[indexMap!(Size, Size, Shape, MatOdr)(i, i)];
	    }
	  }

	  size_t st;
	  foreach(scope i; 0..Size){
	    st= i*Size;
	    result[st..st+Size]= inv2d[i];
	  }
	}
	else{
	  throw new SingularMatrix();
	}

	return result;
      }

      /**********************
       * Lower triangular matrix $(I L)
       **********************/
      T[arrayLength!(Size, Size, MatrixType.lowerTri)]
	matrix(DecomposedMat Mat: DecomposedMat.lowerTri)() nothrow @nogc{
	enum TypeOfIndex LEN= arrayLength!(Size, Size, MatrixType.lowerTri);
	T[LEN] result= void;

	switch(Method){
	case DecompScheme.partialPivDoolittle:
	  {	// diagonal elements
	    auto idxSetDest= IndexSetDiag!(Size, Size, MatrixType.lowerTri, MatOdr)();
	    foreach(scope TypeOfIndex idx; idxSetDest) result[idx]= VALUE_ONE;
	  }
	  {	// strict lower triangular elements
	    auto idxSetSrc= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
	    auto idxSetDest= IndexSetStrictTriL!(Size, MatrixType.lowerTri, MatOdr)();
	    foreach(scope idxSrc, idxDest; zip(idxSetSrc, idxSetDest)) result[idxDest]= _values[idxSrc];
	  }
	  break;
	case DecompScheme.partialPivCrout:
	  auto idxSetSrc0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
	  auto idxSetSrc1= IndexSetStrictTriL!(Size, MatrixType.dense, MatOdr)();
	  foreach(scope idxSrc, idxDest;
		  zip(merge(idxSetSrc0, idxSetSrc1), iota(LEN))) result[idxDest]= _values[idxSrc];
	  break;
	default:
	  assert(false);
	}
	return result;
      }

      /************************
       * Upper triangular matrix $(I U)
       ************************/
      T[arrayLength!(Size, Size, MatrixType.upperTri)]
	matrix(DecomposedMat Mat: DecomposedMat.upperTri)() nothrow @nogc{
	enum TypeOfIndex LEN= arrayLength!(Size, Size, MatrixType.upperTri);
	T[LEN] result= void;

	switch(Method){
	case DecompScheme.partialPivDoolittle:
	  auto idxSetSrc0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
	  auto idxSetSrc1= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
	  foreach(scope idxSrc, idxDest;
		  zip(merge(idxSetSrc0, idxSetSrc1), iota(LEN))
		  ) result[idxDest]= _values[idxSrc];
	  break;
	case DecompScheme.partialPivCrout:
	  {	// diagonal elements
	    auto idxSetDest= IndexSetDiag!(Size, Size, MatrixType.upperTri, MatOdr)();
	    foreach(scope TypeOfIndex idx; idxSetDest) result[idx]= VALUE_ONE;
	  }
	  {	// strict upper triangular elements
	    auto idxSetSrc= IndexSetStrictTriR!(Size, MatrixType.dense, MatOdr)();
	    auto idxSetDest= IndexSetStrictTriR!(Size, MatrixType.upperTri, MatOdr)();
	    foreach(sceop idxSrc, idxDest; zip(idxSetSrc, idxSetDest)) result[idxDest]= _values[idxSrc];
	  }
	  break;
	default:
	  assert(false);
	}
	return result;
      }
    }

  private:
    T[arrayLength!(Size, Size, MatrixType.dense)] _values;
    bool _invertible= false;
  }
};
