/******************************************************************************
 * Module `numeric.sekitk.matirx'
 *
 * for `qvm'
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 0.0
 ******************************************************************************/
module sekitk.matrix;

import std.algorithm: min;

mixin template MatrixImpl(T, real Threshold){
	import sekitk.base;
	import sekitk.exceptions: SingularMatrix;
	import sekitk.indexset: IndexSetDiag,
		IndexSetSubDiag,
		IndexSetStrictTriR,
		IndexSetStrictTriL,
		IndexSetTranspose;
	import sekitk.lu;
	import sekitk.eigen;

	/************************************************************
	 * class of matrix
	 *
	 * Row -> 行數
	 * Column -> 列數
	 * Shape -> 行列要素の形狀 dense, diagonal, tridiagonal, upper triangular, lower triangular.
	 ************************************************************/
	struct Matrix(TypeOfSize Row, TypeOfSize Column,
							 MatrixType Shape: MatrixType.zero,
								MajorOrder MatOdr= MajorOrder.row)
	if(Row > 0 && Column > 0){
		import std.traits: Unqual;

		private{
			alias TypeOfThis= typeof(this);
		}

		@safe pure nothrow @nogc{
			this(){}
			this(in TypeOfThis other){}
		}

		@safe pure nothrow const{
			// sign
			TypeOfThis opUnary(string Op)() @nogc
			if(isPlusOrMinusSign!Op){
				return new typeof(return)();
			}

			// add & sub
		  Matrix!(Row, Column,
							ShapeR,
							MatOdr) opBinary(string Op,
															 MatrixType ShapeR)(in Matrix!(Row, Column,
																														 ShapeR,
																														 MatOdr) rhs) @nogc
			if(isPlusOrMinusSign!Op){
				return new typeof(return)(rhs);
			}

			/****************************************
			 * O * scalar, O / scalar
			 ****************************************/
		  TypeOfThis opBinary(string Op, TypeR)(in TypeR rhs) @nogc
			if((Op == "*" || Op == "/") && is(TypeR: T)){
				return new typeof(return)();
			}

			// scalar * matrix
			// O(m*n)
			TypeOfThis opBinaryRight(string OP: "*")(in T lhs){
				return new typeof(return)();
			}

			// matrix * vector
		  Matrix!(Row, 1, Shape, MatOdr) opBinary(string OP: "*")(in Vector!Column rhs){
			  return new typeof(return)();
			}

			/****************************************
			 * O * Matrix
			 ****************************************/
			Matrix!(Row, ColumnR,
						  Shape, MatOdr) opBinary(string OP: "*",
																			TypeOfSize ColumnR,
																			MatrixType ShapeR)(in Matrix!(Column, ColumnR,
																																		ShapeR, MatOdr) rhs){
				return new typeof(return)();
			}
		}
	}

	/// ditto
	class Matrix(TypeOfSize Row, TypeOfSize Column,
							 MatrixType Shape, MajorOrder MatOdr= MajorOrder.row)
	if(Shape !is MatrixType.zero && matrixConstraint!(Row, Column, Shape)){
		import std.traits: Unqual, TemplateArgsOf,
			isDynamicArray, isStaticArray, isIntegral;
		import std.range: isOutputRange;
		import std.typecons: Tuple, Ternary;
		import std.format: FormatSpec;
		import std.complex;

		import sekitk.base;
		import sekitk.attribute;

		/******************************************
		 * Alias and manifest constants
		 ******************************************/
	private:
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

		enum T VALUE_ZERO= T(0.0L);
		enum bool isSquare= (Row == Column)? true : false;
		static if(isSquare){
			enum TypeOfSize Size= Row;
			alias LibLU= LUdecomposition!(T, Threshold, Size, Shape,
																		(MatOdr is MajorOrder.row)?
																			LUdAlgorithm.doolittle: LUdAlgorithm.crout,
																		MatOdr);
		}

	public:
		/******************************************
		 * mixin
		 ******************************************/
		static if(Shape is MatrixType.dense){
			import sekitk.dense;
			static if(isSquare){
			  mixin(SQ_DENSE);
			}
			else{
				mixin(RECT_DENSE);
			}
		}
		else static if(Shape is MatrixType.band1){
			import sekitk.band1;
			mixin(BAND1_IMPL);
		}
		else static if(Shape is MatrixType.band3){
			import sekitk.band3;
			mixin Band3!(T, Threshold, Size, MatOdr) MethodsImpl;
		}
		else static if(Shape is MatrixType.upperTri){
			import sekitk.uppertri;
			mixin(UPPER_TRI);
		}
		else static if(Shape is MatrixType.lowerTri){
			import sekitk.lowertri;
			mixin(LOWER_TRI);
		}

	  // Constructors
	  @safe pure nothrow{
			import std.algorithm: copy;

			/// Default constructor
			this() @nogc{}

			/****************************************
			 * Copy constructor (the same shape)
			 * Params:
			 * 	other= another matrix of type "Matrix!(Row, Column, Shape, MatOdr)"
			 ****************************************/
			this(inout scope TypeOfThis other){
				this._values[]= other._values[];
				this._attr= other._attr;
				static if(isSquare){
					if(other.ev !is null) this.ev= new typeof(ev)(other.ev);	// eigen values

					final switch(Shape){	// LU decomposition
					case MatrixType.dense, MatrixType.band3:
						if(other.lup !is null) this.lup= new typeof(this.lup)(other.lup);
						break;
					case MatrixType.zero, MatrixType.band1, MatrixType.upperTri, MatrixType.lowerTri:
						break;
					}
				}
			}

			/****************************************
			 * Initialize with a number of type "T"
			 * Params:
			 * 	num= a number
			 ****************************************/
			this(TypeScalar)(inout scope TypeScalar num) @nogc
		  if(is(TypeScalar: T)){
				_values[]= num;
			}

			/****************************************
			 * Initialize with an array
			 * Params:
			 * 	nums= array of type ArrayType
			 * 	attrGiven= 
			 ****************************************/
			this(ArrayType)(scope inout ArrayType nums,
										  scope inout TypeOfAttr attrGiven= TypeOfAttr.init) @nogc
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
				this(ArrayType)(scope inout T[Column][Row] num,
												scope inout TypeOfAttr attrGiven= TypeOfAttr.init) @nogc
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
				this(scope inout Vector!Row[Column] vecArray,
												 scope inout TypeOfAttr attrGiven= TypeOfAttr.init) @nogc{
					static if(MatOdr is MajorOrder.row){
						size_t idx= void;
						foreach(idx_h; 0..Column){
							idx= idx_h;
							foreach(idx_v; 0u..Row){
								this._values[idx]= vecArray[idx_h]._values[idx_v];
								idx += Column;
							}
						}
					}
					else{
						size_t idx_st= 0u, idx_en= Row;
						foreach(idx_h; 0u..Column){
							this._values[idx_st..idx_en]= vecArray[idx_h]._values[];
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
		@safe pure nothrow const{
			TypeOfThis opUnary(string Op: "+")() @nogc{
				return new typeof(return)(this);
			}

			TypeOfThis opUnary(string Op: "-")() @nogc{
				TypeOfInternalArray result= _values[];
				foreach(ref elm; result) elm= -elm;
				return new typeof(return)(result);
			}

			/****************************************
			 * Matrix ± Matrix
			 ****************************************/
		  TypeOfThis opBinary(string Op)(inout scope TypeOfThis rhs) @nogc
			if(isPlusOrMinusSign!Op){
			  TypeOfInternalArray temp= _values[]; 
				mixin("temp[] " ~Op ~"= rhs._values[];");
				return new typeof(return)(temp);
			}

			/// ditto
		  TypeOfThis opBinary(string Op,
													MatrixType ShapeR: MatrixType.zero
													)(inout scope Matrix!(Row, Column, ShapeR, MatOdr) rhs) @nogc
			if(isPlusOrMinusSign!Op && ShapeR != Shape){
				return new typeof(return)(this);
			}


			/****************************************
			 * Matrix * Scalar, Matrix / Scalar
			 ****************************************/
			// matrix * scalar or matrix / scalar
			// O(m*n)
		  TypeOfThis opBinary(string Op, TypeR)(inout scope TypeR rhs) @nogc
			if((Op == "*" || Op == "/") && is(TypeR: T)){
			  TypeOfInternalArray num= _values[]; 
				mixin("num[] " ~Op ~"= rhs;");
				return new typeof(return)(num);
			}

			/****************************************
			 * Scalar * Matrix
			 ****************************************/
			// O(m*n)
		  TypeOfThis opBinaryRight(string Op: "*", TypeL)(inout scope TypeL lhs) @nogc
			if(is(TypeL: T)){
				return this.opBinary!(Op, TypeL)(lhs);
			}

			/****************************************
			 * Matrix * Vector
			 ****************************************/
			Vector!Row opBinary(string OP: "*")(inout scope Vector!Column rhs) @nogc{
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
															 )(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs) @nogc{
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
			Vector!Row opIndex(IdxType)(inout scope IdxType j)
			if(isIntegral!IdxType)
			in(j > 0u && j <= Column){
				return new typeof(return)(sliceColumn(j-1u));
			}

			/****************************************
			 * element
			 * Throws:
			 *  RangeError
			 ****************************************/
			T opIndex(IdxType)(inout scope IdxType idxRow, inout scope IdxType idxColumn) @nogc
			if(isIntegral!IdxType){
				typeof(return) result;
				const auto pos= MatrixPosition!(Row, Column)(idxRow, idxColumn);
				if(pos.rangeCheck){
				  auto map= indexOfInternalArray!(Row, Column, Shape, MatOdr)(pos);
					result= map.isZero? VALUE_ZERO : _values[map.index];
				}
				return result;
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
			Array2D opCast(Array2D: T[Column][Row])() @nogc
			if(((MatOdr is MajorOrder.row) && is(Array2D: T[Column][Row])) ||
				 ((MatOdr is MajorOrder.column) && is(Array2D: T[Row][Column]))){
				enum ITR_LIMIT_OUTER= (MatOdr is MajorOrder.row)? Row: Column;
				enum ITR_LIMIT_INNER= (MatOdr is MajorOrder.row)? Column: Row;
				typeof(return) num= void;

				foreach(i; 0u..ITR_LIMIT_OUTER)
					foreach(j; 0u..ITR_LIMIT_INNER){{
						num[i][j]= this.opIndex(i+1u, j+1u);
			  }}

				return num;
			}
		}	// end of block `@safe pure nothrow const'

	  // Reserved methods
		@safe const{
			/****************************************
			 * Convert to string
			 ****************************************/
			override string toString(){
				char[] buf;
				buf.reserve(6u*Row*Column+(8u+2u*Row)*char.sizeof);
				auto fmt = FormatSpec!char("%s");
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

			  immutable string LEFT_PAREN= "Matrix[";
			  immutable string TAB= array('\x20'.repeat(LEFT_PAREN.length));
				T num= void;
				Tuple!(TypeOfIndex, "index", bool, "isZero") map;
				w.put(LEFT_PAREN);
				foreach(TypeOfSize i; 1u..Row+1u){
					if(i > 1u) put(w, TAB);
					w.put("[");
					foreach(TypeOfSize j; 1u..Column+1u){
						map= indexOfInternalArray!TemplateArgs(MatrixPosition!(Row, Column)(i, j));
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
					w.put("]");
					if(i < Row) w.put("\n");
				}
				w.put("]");
			}
		}

	  // Other common methods
		@safe pure{
			/****************************************
			 * Transpose
			 ****************************************/
			Matrix!(Column, Row,
							MatrixTypeOfTranspose!Shape,
							MatOdr) transpose() nothrow const @property{
				import std.range: iota, zip;
				enum TypeOfIndex LEN= arrayLength!(TemplateArgs[0..3]);
				TypeOfInternalArray num= void;
				auto idxSet= IndexSetTranspose!TemplateArgs();
				foreach(idxTr, idxSrc; zip(idxSet, iota(LEN))) num[idxTr]= _values[idxSrc];
				return new typeof(return)(num);//, _attr);
			}

			/****************************************
			 * Decomposition
			 ****************************************/
			void decompose(DecompScheme Scheme)()
			if(isSquare || (!isSquare && Scheme is DecompScheme.singularValue)){
				final switch(Scheme){
				case DecompScheme.singularValue:
					assert(false, "FIXME: SV decomposition");
					break;
				case DecompScheme.qr:
					assert(false, "FIXME: QR decomposition");
					break;
				case DecompScheme.lup:
					if(lup is null){
						lup= LibLU.partialPivLU(this._values);
					}
					break;
				case DecompScheme.lupq:
					assert(false, "FIXME: LUPQ decomposition");
					break;
				case DecompScheme.ldu:
					assert(false, "FIXME: LDU decomposition");
					break;
				case DecompScheme.cholesky:
					assert(false, "FIXME: Cholesky decomposition");
				}
			}
		}

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

		// Methods of square matrix
		static if(isSquare){
			/****************************************
			 * get the decomposed matrix
			 ****************************************/
			auto decomposedMatrix(DecompScheme Scheme: DecompScheme.lup,
														DecomposedMat Mat)() @safe pure nothrow const
		  if(Mat is DecomposedMat.lowerTri ||
				 Mat is DecomposedMat.upperTri ||
				 Mat is DecomposedMat.permutationLeft){
				static if(Mat is DecomposedMat.lowerTri){
					enum MatrixType ShapeOfReturn= MatrixType.lowerTri;
				}
				else static if(Mat is DecomposedMat.upperTri){
					enum MatrixType ShapeOfReturn= MatrixType.upperTri;
				}
				else{
					enum MatrixType ShapeOfReturn= MatrixType.dense;
				}
				return new Matrix!(Size, Size, ShapeOfReturn, MatOdr)(lup.matrix!Mat);
			}

			@property{
				/**************************************
				 * Trace
				 **************************************/
				T trace() @safe pure nothrow @nogc const{
					import std.range: dropOne;
					import numeric.pseudocmplx: assumeRealNum;
					typeof(return) result= void;

					if(ev is null){
						result= _values[0];
						static if(Size > 1u){
							auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();
							foreach(TypeOfIndex idx; idxSet.dropOne) result += _values[idx];
						}
					}
					else{
						static if(isComplex!T) result= ev.trace;
						else result= assumeRealNum!(CT, Threshold)(ev.trace);
					}
					return result;
				}

				/**************************************
				 * Determinant
				 **************************************/
				static if(Size < 4){
					T det() @safe pure nothrow @nogc const{
						if(_attr.isSingular) return VALUE_ZERO;
						else return detImpl;
					}
				}
				else{
					T det() @safe pure nothrow{
						import numeric.pseudocmplx: assumeRealNum;
						typeof(return) result= void;

						if(_attr.isSingular){
							result= VALUE_ZERO;
						}
						else{
						  if(ev){
								static if(isComplex!T) result= ev.det;
								else result= assumeRealNum!(CT, Threshold)(ev.det);
							}
							else result= detImpl;	// if DIP 1008 is enable, @nogc 
						}

						return result;
					}
				}

				/**************************************
				 * Inverse matrix
				 **************************************/
				typeof(this) inverse() @safe pure{
				//if(!(flags & CHECK_SINGULARITY)) {}//checkSingularity(lup);	// fwelm or lup or ldu
				//if(flags & TYP_INVERTIBLE){
					return new typeof(return)(inverseImpl);
					/*}
				else{
					throw new SingularMatrix;
					}*/
				}
			}

			/****************************************
			 * Eigen value
			 *
			 * \Lambda= \lambda_1, \lambda_2, ..., \lambda_n
			 ****************************************/
			void findEigenValues() @safe @property{//@safe pure nothrow @property{
				import std.stdio;	// DEBUG:
				if(ev is null){
					final switch(Shape){
					case MatrixType.zero:
						ev= new typeof(ev);
						break;
					case MatrixType.dense:
						static if(Size == 1u){
							ev= new typeof(ev)(_values);
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
							import mathematic.equations: QuadraticEq;
							auto eq= new QuadraticEq!(T, Threshold)(1.0, -trace, det);
							ev= new typeof(ev)(eq.solve);
						}
						else{	// n > 2
							//if(qr is null) decompose(DecompMethod.QR);
							assert(false, "QR decomposition is not available.");	// FIXME: DO SOMETHING
						}
						break;
					case MatrixType.band3:
						assert(false, "FIXME: eigen value computation of a tridiagonal matrix");
						break;
					case MatrixType.band1:
						static if(is(T == RT)){
							CT[Size] temp= void;
							foreach(idx; 0u..Size) temp[idx]= complex!RT(_values[idx]);
							ev= new typeof(ev)(temp);
						}
						else ev= new typeof(ev)(_values);
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
						ev= new typeof(ev)(num);
					}
				}
				else{}	// NOP: The eigen values have been already obtained
			}

			/****************************************
			 * eigen values
			 ****************************************/
			CT[Size] eigenValues() @safe {//@safe pure nothrow @property{
				if(ev is null) findEigenValues;
				return ev.values;
			}
		}	// end of block "static if(isSquare)"

	package:
		TypeOfInternalArray _values;
		TypeOfAttr _attr;
		static if(isSquare){
			enum LUdAlgorithm AlgoLU= LUdAlgorithm.doolittle;
			Eigen!(CT, Threshold, Size) ev;
		}

	public:
		static if(isSquare) LibLU.ArrayLUP lup;

  private:
		pure const{
			bool approxEqualZero(inout scope T num) @safe nothrow @nogc{
				import numeric.approx: approxEqualZero;
				return approxEqualZero!(T, Threshold)(num);
			}

			auto trustedAssumeUnique(U)(U t) @trusted nothrow @nogc{
				import std.exception: assumeUnique;
				return assumeUnique(t);
			}

		// returns vector
			T[Row] sliceColumn(IdxType)(in IdxType col) @safe nothrow @nogc
			if(isIntegral!IdxType)
		  in(col >= 0u && col < Column){
			  typeof(return) temp= void;

				final switch(Shape){
				case MatrixType.zero:
					temp[]= VALUE_ZERO;
					break;
				case MatrixType.dense:
					final switch(MatOdr){
					case MajorOrder.column:
						size_t st= Row*col;
						temp[]= _aluesv[st .. st+Row];
						break;
					case MajorOrder.row:
						foreach(i; 0u..Row) temp[i]= this.opIndex(i+1u, col+1u);
					}
					break;
				case MatrixType.band1:
					foreach(i; 0u..Row){
						if(i == col) temp[i]= _values[i];
						else temp[i]= VALUE_ZERO;
					}
					break;
				case MatrixType.band3, MatrixType.upperTri, MatrixType.lowerTri:
					foreach(i; 0u..Row) temp[i]= this.opIndex(i+1u, col+1u);
				}
				return temp;
		  }
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
