/*****************************************************************************
 * This module is for internal use only.
 *
 * The implementations of all type of matrix are combined to "Matrix" here.
 * If you use this directly, it would not be useful.
 *
 * Version: 1.0
 * License: Boost License 1.0
 * Authors: 新井 浩平 (Kohei ARAI), arawi_kohei_takasaki@yahoo.co.jp
 * Source: https://github.com/ARAWI-Kohei-tgj/sekitk/src/qvm/matrix.d
 *****************************************************************************/
module sekitk.qvm.matrix;

import std.traits: isFloatingPoint;
import sekitk.complex.pseudo: isComplex;

mixin template MatrixImpl(T, real Threshold) if(isFloatingPoint!T || isComplex!T){
  import sekitk.qvm.common: TypeOfSize, TypeOfIndex, matrixConstraint;
  import sekitk.qvm.impl.zero, sekitk.qvm.impl.permutation;

  /***********************************************************
   * matrix
   *
   * This type can be customized by setting its size, matrix type and major order
   * via the template arguments.
   *
   * Template arguments:
   * Shape
   * MatrixType.zero= zero matrix
   * MatrixType.permutation= permutation matrix
   * MatrixType.dense= dense matrix
   * MatrixType.band1= diagonal matrix
   * MatrixType.band3= tridiagonal matrix
   * MatrixType.upperTri= upper triangular matrix
   * MatrixType.lowerTri= lower triangular matrix
   *
   * MatOdr
   * MajorOrder.row= row major order
   * MajorOrder.column= column major order
   ***********************************************************/
  mixin(MATRIX_IMPL_ZERO);
  mixin(MATRIX_IMPL_PERM);

  // dense, diagonal, tridiagonal, upper triangular, lower triangular
  class Matrix(TypeOfSize Row, TypeOfSize Column,
	       MatrixType Shape, MajorOrder MatOdr= MajorOrder.row)
  if(Shape !is MatrixType.zero &&
     Shape !is MatrixType.permutation &&
     matrixConstraint!(Row, Column, Shape)){
    import std.complex;
    import std.format: FormatSpec;
    import std.meta: Alias, AliasSeq;
    import std.range: isOutputRange;
    import std.traits: Unqual, TemplateArgsOf,
      isDynamicArray, isStaticArray, isIntegral;
    import std.typecons: Tuple, Ternary, Nullable;
    import std.variant: VariantN, maxSize;

    import sekitk.qvm.common;
    import sekitk.qvm.exception: SingularMatrix;
    import sekitk.qvm.indexset: IndexSetDiag,
      IndexSetSubDiag,
      IndexSetStrictTriR,
      IndexSetStrictTriL,
      IndexSetTranspose;
    import sekitk.qvm.algorithm.decomposition: MatrixDecomp;
    import sekitk.qvm.eigen;
    import sekitk.qvm.attribute;

  private:
    /*****************************************
     * Alias and manifest constants
     *****************************************/
    static if(isFloatingPoint!T){
      alias RT= T;	// type of real number
      alias CT= Complex!T;	// type of complex number
    }
    else{
      alias RT= BaseRealType!T;
      alias CT= T;
    }

    alias TypeOfThis= typeof(this);
    alias TemplateArgs= TemplateArgsOf!(typeof(this));
    alias TypeOfInternalArray= T[arrayLength!(TemplateArgs[0..3])];
    alias TypeOfAttr= MatrixAttribute!(TemplateArgs[0..3]);

    enum T VALUE_ZERO= Identity!(T, "+");
    enum bool isSquare= (Row == Column);
    static if(isSquare){
      enum TypeOfSize Size= Row;
      //alias LibLU= MatrixDecomp!(T, Threshold).DecompLU!(Size, Shape, MatOdr);
    }

    public @safe:
    /******************************************
     * mixin
     ******************************************/
    static if(Shape is MatrixType.dense){
      import sekitk.qvm.impl.dense;
      static if(isSquare){
	mixin(SQ_DENSE);
      }
      else{
	mixin(RECT_DENSE);
      }
    }
    else static if(Shape is MatrixType.band1){
      import sekitk.qvm.impl.band1;
      mixin(BAND1_IMPL);
    }
    else static if(Shape is MatrixType.band3){
      import sekitk.qvm.impl.band3;
      mixin Band3!(T, Threshold, Size, MatOdr) MethodsImpl;
    }
    else static if(Shape is MatrixType.upperTri){
      import sekitk.qvm.impl.uppertri;
      mixin(UPPER_TRI);
    }
    else static if(Shape is MatrixType.lowerTri){
      import sekitk.qvm.impl.lowertri;
      mixin(LOWER_TRI);
    }
    else{
      static assert(false);
    }

    // Constructors
    pure nothrow @nogc{
      /// Default constructor
      this(){}

      /****************************************
       * Copy constructor (the same shape)
       * Params:
       * 	other= another matrix of type "Matrix!(Row, Column, Shape, MatOdr)"
       ****************************************/
      this(ref return scope inout TypeOfThis other){
	_values[]= other._values[];
	_attr= other._attr;
      }

      /****************************************
       * Initialize with a number of type "T"
       * Params:
       * 	num= a number
       ****************************************/
      this(TypeScalar)(scope const TypeScalar num)
      if(is(TypeScalar: T)){
	_values[]= num;
      }

      /****************************************
       * Initialize with an array
       * Params:
       * 	nums= array of type ArrayType
       * 	attrGiven= 
       ****************************************/
      this(ArrayType)(scope const ArrayType nums,
		      scope const TypeOfAttr attrGiven= TypeOfAttr.init)
      if((isDynamicArray!ArrayType && is(Unqual!ArrayType: T[])) ||
	 (isStaticArray!ArrayType && is(Unqual!ArrayType: TypeOfInternalArray)))
      in(nums.length == arrayLength!(TemplateArgs[0..3]), "mismatch array length"){
	  _values[]= nums[];
	  _attr= attrGiven;
      }

      // for only dense matrix
      static if(Shape is MatrixType.dense){
	/**************************************
	 * Initialize with a 2-dimensional array
	 **************************************/
	this(ArrayType)(scope const T[Column][Row] num,
			scope const TypeOfAttr attrGiven= TypeOfAttr.init)
	if((MatOdr is MajorOrder.row) && is(Unqual!ArrayType: T[Column][Row]) ||
	   (MatOdr is MajorOrder.column) && is(Unqual!ArrayType: T[Row][Column])){
	  enum TypeOfSize ITR_LIMIT= (MatOdr is MajorOrder.row)? Row: Column;
	  enum TypeOfSize STEP= (MatOdr is MajorOrder.row)? Column: Row;
	  size_t st, en= 0u;

	  foreach(i; 0u..ITR_LIMIT){
	    st= en;
	    en += STEP;
	    _values[st..en]= num[i][];
	  }
	  _attr= attrGiven;
	}

	/**************************************
	 * Initialize with an array of Vector!Row
	 **************************************/
	this(scope const Vector!Row[Column] vecArray,
	     scope const TypeOfAttr attrGiven= TypeOfAttr.init){
	  static if(MatOdr is MajorOrder.row){
	    size_t idx= void;
	    foreach(idxH; 0u..Column){
	      idx= idxH;
	      foreach(idxV; 0u..Row){
		this._values[idx]= vecArray[idxH]._values[idxV];
		idx += Column;
	      }
	    }
	  }
	  else{
	    size_t idx_st= 0u, idx_en= Row;
	    foreach(idxH; 0u..Column){
	      this._values[idx_st..idx_en]= vecArray[idxH]._values[];
	      idx_st= idx_en;
	      idx_en += Row;
	    }
	  }
	  _attr= attrGiven;
	}
	
      }	// end of block "Shape is MatrixType.dense"
    }

    // Operators & Reserved methods
    /// plus or minus sign
    pure nothrow const{
      TypeOfThis opUnary(string Op: "+")(){
	return new typeof(return)(this);
      }

      TypeOfThis opUnary(string Op: "-")(){
	TypeOfInternalArray result= _values[];
	foreach(ref elm; result) elm= -elm;
	return new typeof(return)(result);
      }

      /****************************************
       * Matrix ± Matrix
       ****************************************/
      TypeOfThis opBinary(string Op)(inout scope TypeOfThis rhs)
	if(isPlusOrMinusSign!Op){
	  TypeOfInternalArray temp= _values[]; 
	  mixin("temp[] " ~Op ~"= rhs._values[];");
	  return new typeof(return)(temp);
	}

      /// ditto
      TypeOfThis opBinary(string Op,
			  MatrixType ShapeR: MatrixType.zero
			  )(inout scope Matrix!(Row, Column, ShapeR, MatOdr) rhs)
      if(isPlusOrMinusSign!Op && ShapeR != Shape){
	return new typeof(return)(this);
      }

      /****************************************
       * Matrix * Scalar, Matrix / Scalar
       ****************************************/
      // matrix * scalar or matrix / scalar
      // O(m*n)
      TypeOfThis opBinary(string Op, TypeR)(inout scope TypeR rhs)
      if((Op == "*" || Op == "/") && is(TypeR: T)){
	TypeOfInternalArray num= _values[]; 
	mixin("num[] " ~Op ~"= rhs;");
	return new typeof(return)(num);
      }

      /****************************************
       * Scalar * Matrix
       ****************************************/
      // O(m*n)
      TypeOfThis opBinaryRight(string Op: "*", TypeL)(inout scope TypeL lhs)
      if(is(TypeL: T)){
	return this.opBinary!(Op, TypeL)(lhs);
      }

      /****************************************
       * Matrix * Vector
       ****************************************/
      Vector!Row opBinary(string OP: "*")(inout scope Vector!Column rhs){
	auto colmat= Matrix!(Column, 1u, MatrixType.dense, MatOdr)(rhs.a);
	return new typeof(return)(this*colmat);
      }

      /****************************************
       * Matrix * O
       ****************************************/
      Matrix!(Row, ColumnR,
	      MatrixType.zero,
	      MatOdr) opBinary(string Op: "*",
			       TypeOfSize ColumnR,
			       MatrixType ShapeR: MatrixType.zero
			       )(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs){
	return new typeof(return)();
      }

      /****************************************
       * Index operator
       *
       * Params:
       *  j= index of the column, (1 ≤ j ≤ COLUMN)
       *
       * Returns:
       *  Sliced column as a Vector!ROW
       *
       * Throws:
       *  RangeError
       ****************************************/
      Vector!Row opIndex(IdxType)(in IdxType j) @nogc
      if(isIntegral!IdxType)
      in(j > 0 && j <= Column){
	return typeof(return)(sliceVertically(j-1u));
      }

      /****************************************
       * element
       * Throws:
       *  RangeError
       ****************************************/
      T opIndex(IdxType)(in IdxType idxRow, in IdxType idxColumn) @nogc
      if(isIntegral!IdxType)
      in(idxRow > 0 && idxRow <= Row
	 && idxColumn > 0 && idxColumn <= Column){
	const pos= MatrixPosition!(Row, Column)(idxRow, idxColumn);
        pos.rangeCheck;
	pragma(msg, "through 0");
	auto map= internalIndexOf!(Row, Column, Shape, MatOdr)(pos);
	pragma(msg, "through 1");
	return map.isZero? VALUE_ZERO: _values[map.index];
      }

      /****************************************
       * Cast to dense matrix
       ****************************************/
      Mat opCast(Mat: Matrix!(Row, Column, MatrixType.dense, MatOdr))(){
	typeof(return) result;
	static if(Shape is MatrixType.dense){
	  result= new typeof(return)(this);
	}
	else{
	  result= new typeof(return)(arrayDense);
	}

	return result;
      }

      /****************************************
       * Cast to 2-dimensional array
       ****************************************/
      Array2D opCast(Array2D)() @nogc
      if(is(Array2D: TypeOfInternalArray)){
	typeof(return) num= void;

	foreach(i; 0u..Row)
	  foreach(j; 0u..Column){{
	      num[i][j]= this.opIndex(i+1u, j+1u);
        }}

	return num;
      }
    }	// end of block `@safe pure nothrow const'

    // Reserved methods
    const{
      /****************************************
       * Convert to string
       ****************************************/
      override string toString(){
	string trustedAssumeUnique(char[] bufMutable) @trusted pure{
	  import std.exception: assumeUnique;
	  return assumeUnique(bufMutable);
	}

	char[] buf;
	{
	  enum size_t RESERVE_SIZE= (LEFT_PAREN.length	// "Matrix["
				     +"[],\n".length	// "[],\n" or "[]]\n"
				     +DELIM.length*(Column-1)	// ", "
				     +"-123.45678".sizeof*Column
				     )*Row;
	  buf.reserve(RESERVE_SIZE);
	}
	auto fmt= FormatSpec!char("%s");
	this.toString((const(char)[] s){buf ~= s;}, fmt);

	return trustedAssumeUnique(buf);
      }

      /// ditto
      void toString(Writer, Char)(scope Writer w, FormatSpec!Char formatSpec)
	if(isOutputRange!(Writer, const(Char)[])){
	  import std.array: array;
	  import std.format: formatValue;
	  import std.math: signbit;
	  import std.range: put, repeat;
	  import std.typecons: Tuple;

	  immutable string TAB= array('\x20'.repeat(LEFT_PAREN.length));
	  T num= void;
	  Tuple!(TypeOfIndex, "index", bool, "isZero") map;

	  w.put(LEFT_PAREN);
	  foreach(TypeOfSize i; 1u..Row+1u){
	    if(i > 1) w.put(TAB);
	    w.put("[");
	    foreach(TypeOfSize j; 1u..Column+1u){
	      map= internalIndexOf!TemplateArgs(MatrixPosition!(Row, Column)(i, j));	// FIXME:
	      if(map.isZero){
		w.put("*");	// always zero
	      }
	      else{
		num= _values[map.index];
		if(approxEqualZero(num)) w.put("≈0");
		else w.formatValue(num, formatSpec);
	      }
	      if(j < Column) w.put(", ");
	    }
	    w.put((i < Row)? "],\n": "]");
	  }
	  w.put("]");
	}
    }

    // Other mathematical methods
    pure{
      /****************************************
       * Transpose
       *
       * \bvec{A}^T
       ****************************************/
      Matrix!(Column, Row,
	      ReturnTypeOfTranspose!Shape,
	      MatOdr) transpose() nothrow const @property{
	import std.range: iota, zip;
	enum TypeOfIndex LEN= arrayLength!(TemplateArgs[0..3]);

	TypeOfInternalArray num= void;
	auto idxSet= IndexSetTranspose!TemplateArgs();
	foreach(idxTr, idxSrc; zip(idxSet, iota(LEN))) num[idxTr]= _values[idxSrc];
	return new typeof(return)(num, _attr.transpose);
      }

      /****************************************
       * Adjoint matrix
       *
       *
       ****************************************/
      Matrix!(Column, Row,
	      ReturnTypeOfTranspose!Shape,
	      MatOdr) adjoint() nothrow const @property{
	import std.range: iota, zip;
	static if(isFloatingPoint!T) import sekitk.complex.pseudo: conj;
	enum TypeOfIndex LEN= arrayLength!(TemplateArgs[0..3]);

	TypeOfInternalArray num= void;
	auto idxSet= IndexSetTranspose!TemplateArgs();
	foreach(idxTr, idxSrc; zip(idxSet, iota(LEN))) num[idxTr]= _values[idxSrc].conj;
	return new typeof(return)(num, _attr.transpose);
      }


      /****************************************
       * Adjugate matrix
       *
       * \adj_ij(\bvec{A})
       ****************************************/
/+
      Matrix!(Row-1, Column-1,
              MatrixType.dense, MatOdr) adj(IdxType)(in IdxType idxRow, in IdxType idxColumn)
      if(isIntegral!IdxType)
      in(idxRow > 0 && idxRow <= Row)
      in(idxColumn > 0 && idxColumn <= Column){
        import std.range: iota, chain;

	T[arrayLength!(Row-1, Column-1, MatrixType.dense)] num= void;
	auto num2Dim= opCast!(T[Column][Row])();
	TypeOfIndex k= 0;
	foreach(i; chain(iota(idxRow-1), iota(idxRow, Row))){
	  foreach(j; chain(iota(idxColumn), iota(idxColumn, Column))) num[k++]= num2Dim[i][j];
	}
	return new typeof(return)(num);
      }
+/
    }

    // non-mathematical methods
    pure{
      /****************************************
       * Decomposition
       ****************************************/
      void decompose(DecompScheme Scheme)()
	if(isSquare || (!isSquare && Scheme is DecompScheme.singularValue)){
/+
        static if(isLU!Scheme){
          resultLU.get= LibLU.decompose!Scheme(cast(T[Column][Row])this);
	}
	else{
	assert(false);
	}
	+/
	}
    }
/+
    /******************************************
     * Decomposed matrix
     ******************************************/
    auto decomposedMatrix(DecompScheme Scheme: DecompScheme.singularValue, DecomposedMat Mat)()
    if(Mat is DecomposedMat.unitaryLeft ||
       Mat is DecomposedMat.diagonal ||
       Mat is DecomposedMat.unitaryRight){
      static if(Mat is DecomposedMat.unitaryLeft){
        auto MatTmplParams= Tuple!(Row, Row, MatrixType.dense, MatOdr);
	T[arrayLength!(MatTmplParams[0..3])] num;
	return new Matrix!MatTmplParams(num);
      }
      else static if(Mat is DecomposedMat.diagonal){
        auto MatTmplParams= Tuple!(Row, Column, MatrixType.band1, MatOdr);
        T[arrayLength!(MatTmplParams[0..3])] num;
        return new Matrix!(MatTmplParams)(num);
      }
      else{
        auto MatTmplParams= Tuple!(Column, Column, MatrixType.dense, MatOdr);
        T[arrayLength!(MatTmplParams[0..3])] num;
        return new Matrix!(MatTmplParams)(num);
      }
    }

    auto decomposedMatrix(DecompScheme Scheme: DecompScheme.qr, DecomposedMat Mat)()
    if(Mat is DecomposedMat.unitary ||
       Mat is DecomposedMat.upperTri){
      static if(Mat is DecomposedMat.unitary){
        return new Matrix!(Row, Column, MatrixType.dense, MatOdr)(num);
      }
      else{
	return new Matrix!(Row, Column, MatrixType.upperTri)(num);
      }
    }
+/
    // Methods of square matrix
    static if(isSquare){
/+
      /****************************************
       * get the decomposed matrix
       ****************************************/
      auto decomposedMatrix(DecompScheme Scheme: DecompScheme.lup,
                            DecomposedMat Mat)() @safe pure nothrow const
      if(Mat is DecomposedMat.lowerTri ||
         Mat is DecomposedMat.upperTri ||
	 Mat is DecomposedMat.permutationLeft){
	static if(Mat is DecomposedMat.lowerTri){
					static if(Shape is MatrixType.dense || Shape is MatrixType.band3){
						enum MatrixType ShapeDest= MatrixType.lowerTri;
					}
					else{
						enum MatrixType ShapeDest= MatrixType.band1;
					}
				}
				else static if(Mat is DecomposedMat.upperTri){
					static if(Shape is MatrixType.dense || Shape is MatrixType.band3){
						enum MatrixType ShapeDest= MatrixType.upperTri;
					}
					else{
						enum MatrixType ShapeDest= MatrixType.band1;
					}
				}
				else{
					enum MatrixType ShapeDest= MatrixType.permutation;
				}

				return new Matrix!(Size, Size, ShapeDest, MatOdr)(lup.matrix!Mat);
			}
+/
		  @property{
		    /**************************************
		     * Trace
		     *
		     * \tr(\bvec{A})
		     **************************************/
		    T trace() pure nothrow @nogc const{
		      import std.range: dropOne;
		      typeof(return) result= void;

		      result= _values[0];
		      static if(Size > 1u){
			auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();
			foreach(TypeOfIndex idx; idxSet.dropOne) result += _values[idx];
		      }
		      return result;
		    }

		    /**************************************
		     * Determinant
		     **************************************/
		    T det() pure nothrow @nogc const{
		      return _attr.isSingular? VALUE_ZERO: detImpl;
		    }

		    /**************************************
		     * Inverse matrix
		     **************************************/
		    typeof(this) inverse() pure{
		      if(_attr.isSingular){
			throw new SingularMatrix();
		      }
		      else{
			return new typeof(return)(inverseImpl);
		      }
		    }
		  }

		  /****************************************
		   * Eigen value
		   *
		   * \Lambda= \lambda_1, \lambda_2, ..., \lambda_n
		   ****************************************/
			/+
			void findEigenValues() @property{//@safe pure nothrow @property{
				import std.stdio;	// DEBUG:
				if(!_ev.isNull){
					final switch(Shape){
					case MatrixType.dense:
						static if(Size == 1u){
							_ev.get= Eigen(CT, Threshold, Size)(_values);
						}
						else static if(Size == 2u){
							/****************
							 * 2×2 matrix
							 *
							 * \left\{
							 * \begin{array}
							 * 	trace(M)= \sum_J \lambda_i
							 * 	det(M)= \pi_J \lambda_i
							 * \end{array}
							 * \right.
							 *
							 * index J= {j \in \mathbf{Z}| j= 1,~2,~\cdots,~n}
							 *
							 * trace(M)= \lambda_1 + \lambda_2
							 * det(M)= \lambda_1 \lambda_2
							 * \Lambda= ax^2+bx+c
							 * \Lambda= \dfrac{-b \pm \sqrt{b^2-4ac}}{2*a}
							 * a= 1, b= -trace(M), c= det(M)
							 ****************/
							import sekitk.equations: QuadraticEq;
							auto eq= new QuadraticEq!(T, Threshold)(1.0, -trace, det);
							_ev.get= new Eigen!(CT, Threshold, Size)(eq.solve);
						}
						else{	// n > 2
							//if(qr is null) decompose(DecompMethod.QR);
							assert(false, "QR decomposition is not available.");	// FIXME: DO SOMETHING
						}
						break;
					case MatrixType.band1:
						static if(is(T == RT)){
							CT[Size] temp= void;
							foreach(idx; 0u..Size) temp[idx]= complex!RT(_values[idx]);
							_ev.get= Eigen!(CT, Threshold, Size)(temp);
						}
						else{
							_ev.get= Eigen!(CT, Threshold, Size)(_values);
						}
						break;
					case MatrixType.band3:
						assert(false, "FIXME: eigen value computation of a tridiagonal matrix");
						break;
					case MatrixType.upperTri, MatrixType.lowerTri:
						CT[Size] num= void;
						auto idxSet= IndexSetDiag!TemplateArgs();
						static if(is(T == RT)){
							foreach(TypeOfIndex idx; idxSet) num[idx]= complex!RT(_values[idx]);
						}
						else{
							foreach(TypeOfIndex idx; idxSet) num[idx]= _values[idx];
						}
						_ev.get= Eigen!(CT, Threshold, Size)(num);
						break;
					case MatrixType.zero, MatrixType.permutation:
						assert(false);
					}
				}
				else{}	// NOP: The eigen values have been already obtained
			}
+/
			/****************************************
			 * eigen values
			 ****************************************/
/+
			CT[Size] eigenValues() {//@safe pure nothrow @property{
				if(!_ev.isNull) findEigenValues;
				return _ev.get.values;
			}
+/
		}	// end of block "static if(isSquare)"

  package:
    TypeOfInternalArray _values;
    TypeOfAttr _attr;
		/+
		static if(isSquare){
			alias AllTypesOfResultLU= AliasSeq!(LibLU.Result!(DecompScheme.luWithNonPivDoolittle),
																					//LibLU.Result!(DecompScheme.luWithNonPivCrout),
																					LibLU.Result!(DecompScheme.luWithPartialPivDoolittle),
																					LibLU.Result!(DecompScheme.luWithPartialPivCrout)//,
																					//LibLU.Result!(DecompScheme.luWithFullPivDoolittle),
																					//LibLU.Result!(DecompScheme.luWithFullPivCrout)
																					);
			alias ResultsLU= VariantN!(maxSize!AllTypesOfResultLU, AllTypesOfResultLU);
			Nullable!ResultsLU resultLU;
			Nullable!(Eigen!(CT, Threshold, Size)) _ev;
		}
+/

    private @safe pure nothrow @nogc const:
    enum string LEFT_PAREN= "Matrix[";
    enum string DELIM= ", ";

    bool approxEqualZero(in T num){
      import sekitk.approx: approxEqualZero;
      return approxEqualZero!(T, Threshold)(num);
    }

    package @safe pure nothrow @nogc const:
    /*****************************************
     * extracts the specified row
     *
     * Params:
     *	idxRow= 0-based index of row (0 ≤ idxRow < Row)
     *****************************************/
    T[Column] sliceHorizontally(IdxType)(in IdxType idxRow)
    if(isIntegral!IdxType)
    in(idxRow >= 0 && idxRow < Row, "Argument is out of range."){
      typeof(return) result= void;

      final switch(MatOdr){
      case MajorOrder.row:
	enum size_t IDX_ST= Column*idxRow;
	result[]= _values[IDX_ST .. IDX_ST+Row];
	break;
      case MajorOrder.column:
	foreach(i; 0u..Column) result[i]= this.opIndex(i+1u, idxRow+1u);
      }
      return result;
    }

    /*****************************************
     * extracts the specified column
     *
     * Params:
     * 	idxColumn= 0-based index of column (0 ≤ idxColumn < Column)
     *****************************************/
    T[Row] sliceVertically(IdxType)(in IdxType idxColumn)
    if(isIntegral!IdxType)
    in(idxColumn >= 0 && idxColumn < Column, "Argument is out of range."){
      typeof(return) result= void;

      final switch(MatOdr){
      case MajorOrder.row:
	foreach(i; 0u..Row) result[i]= this.opIndex(i+1u, idxColumn+1u);
	break;
      case MajorOrder.column:
	const size_t IDX_ST= Row*idxColumn;
	result[]= _values[IDX_ST .. IDX_ST+Row];
      }
      return result;
    }

  } // end of class "Matrix!(size, size, shape)"

  /************************************************************
   * tie four matrix to one
   ************************************************************/
  Matrix!(RowUp+RowLw, ColumnL+ColumnR,
	  MatrixType.dense,
	  MatOdr) tie(TypeOfSize RowUp, TypeOfSize RowLw,
		      TypeOfSize ColumnL, TypeOfSize ColumnR,
		      MatrixType ShapeLT, MatrixType ShapeRT,
		      MatrixType ShapeLB, MatrixType ShapeRB,
		      MajorOrder MatOdr
		      )(in Matrix!(RowUp, ColumnL, ShapeLT, MatOdr) matLT,
			in Matrix!(RowUp, ColumnR, ShapeRT, MatOdr) matRT,
			in Matrix!(RowLw, ColumnL, ShapeLB, MatOdr) matLB,
			in Matrix!(RowLw, ColumnR, ShapeRB, MatOdr) matRB
			) @safe pure nothrow @nogc{
    enum TypeOfSize Row= RowUp+RowLw;
    enum TypeOfSize Column= ColumnL+ColumnR;
    enum MatrixType ShapeOfReturn= (){
      MatrixType result;
      static if(ShapeRT is MatrixType.zero && ShapeLB is MatrixType.zero){
	static if(ShapeLT is ShapeRB &&
		  (ShapeLT is MatrixType.band1 || ShapeLT is MatrixType.band3)){
	  result= ShapeLT;
	}
	else static if((ShapeLT is MatrixType.band1 && ShapeRB is MatrixType.band3) ||
		       (ShapeLT is MatrixType.band3 && ShapeRB is MatrixType.band1)){
	  result= MatrixType.band3;
	}
      }
      else static if(ShapeLT is MatrixType.upperTri &&
		     ShapeLB is MatrixType.zero &&
		     ShapeRB is MatrixType.upperTri){
	result= MatrixType.upperTri;
      }
      else static if(ShapeLT is MatrixType.lowerTri &&
		     ShapeRT is MatrixType.zero &&
		     ShapeRB is MatrixType.lowewrTri){
	result= MatrixType.lowerTri;
      }
      else{
	result= MatrixType.dense;
      }
      return result;
    }();

    const T[ColumnL][RowUp] aryLT= cast(T[ColumnL][RowUp])matLT;
    const T[ColumnR][RowUp] aryRT= cast(T[ColumnR][RowUp])matRT;
    const T[ColumnL][RowLw] aryLB= cast(T[ColumnL][RowLw])matLB;
    const T[ColumnR][RowLw] aryRB= cast(T[ColumnR][RowLw])matRB;
    T[arrayLength!(TemplateArgsOf!(typeof(return))[0..3])] num;

    foreach(i; 0u..RowUp){
      foreach(j; 0u..ColumnL) num[Column*i+j]= aryLT[i][j];
      foreach(j; 0u..ColumnR) num[Column*i+ColumnL+j]= aryRT[i][j];
    }

    foreach(i; 0u..RowLw){
      foreach(j; 0u..ColumnL) num[Column*(RowUp+i)+j]= aryLB[i][j];
      foreach(j; 0u..ColumnR) num[Column*(RowUp+i)+ColumnL+j]= aryRB[i][j];
    }

    return new typeof(return)(num);
  }
}
