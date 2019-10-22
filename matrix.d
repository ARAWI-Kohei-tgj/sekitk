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
				return typeof(return)();
			}

			// add & sub
		  Matrix!(Row, Column,
							ShapeR,
							MatOdr) opBinary(string Op,
															 MatrixType ShapeR)(in Matrix!(Row, Column,
																														 ShapeR,
																														 MatOdr) rhs) @nogc
			if(isPlusOrMinusSign!Op){
				return typeof(return)(rhs);
			}

		  TypeOfThis opBinary(string Op)(in T rhs) @nogc
			if(Op == "*" || Op == "/"){
				return typeof(return)();
			}

			// scalar * matrix
			// O(m*n)
			TypeOfThis opBinaryRight(string OP: "*")(in T lhs){
				return typeof(return)();
			}

			// matrix * vector
		  Matrix!(Row, 1, Shape, MatOdr) opBinary(string OP: "*")(in Vector!Column rhs){
			  return typeof(return)();
			}

			// matrix multiplication
			Matrix!(Row, ColumnR,
						  Shape, MatOdr) opBinary(string OP: "*",
																			TypeOfSize ColumnR,
																			MatrixType ShapeR)(in Matrix!(Column, ColumnR,
																																		ShapeR, MatOdr) rhs){
				return typeof(return)();
			}
		}
	}

	/// ditto
	class Matrix(TypeOfSize Row, TypeOfSize Column,
							 MatrixType Shape, MajorOrder MatOdr= MajorOrder.row)
	if(Shape !is MatrixType.zero && matrixConstraint!(Row, Column, Shape)){
		import std.traits: Unqual, TemplateArgsOf;
		import std.range: isOutputRange;
		import std.typecons: Tuple, Ternary;
		import std.format: FormatSpec;
		import std.complex;
		import sekitk.base;

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
		alias TmplArgsOfThis= TemplateArgsOf!(typeof(this));
		alias TypeOfInternalArray= T[arrayLength!(Row, Column, Shape)];

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
				mixin SqDense!(T, Threshold, Size, MatOdr) MethodsImpl;
			}
			else{
				mixin RectDense!(T, Threshold, Row, Column, MatOdr) MethodsImpl;
			}
		}
		else static if(Shape is MatrixType.band1){
			import sekitk.band1;
			mixin Band1!(T, Threshold, Size, MatOdr) MethodsImpl;
		}
		else static if(Shape is MatrixType.band3){
			import sekitk.band3;
			mixin Band3!(T, Threshold, Size, MatOdr) MethodsImpl;
		}
		else static if(Shape is MatrixType.upperTri){
			import sekitk.uppertri: UpperTriangular;
			mixin UpperTriangular!(T, Threshold, Size, MatOdr) MethodsImpl;
		}
		else static if(Shape is MatrixType.lowerTri){
			import sekitk.lowertri: LowerTriangular;
			mixin LowerTriangular!(T, Threshold, Size, MatOdr) MethodsImpl;
		}

		/******************************************
		 * Constructors
		 ******************************************/
	  @safe pure nothrow{
			import std.algorithm: copy;

			/// Default constructor
			this() @nogc{
			  attrInitialize();
			}

			/************************
			 * Copy constructor (the same shape)
			 * Params:
			 * 	other= another matrix of type "Matrix!(Row, Column, Shape, MatOdr)"
			 ************************/
			this(in TypeOfThis other) @nogc{
				this.v[]= other.v[];

				static if(isSquare){
					this.attr.isInvertible= other.attr.isInvertible;
					if(other.ev !is null) this.ev= new typeof(ev)(other.ev);	// eigen values

					final switch(Shape){	// LU decomposition
					case MatrixType.dense, MatrixType.band3:
						if(other.lup !is null) this.lup= new typeof(this.lup)(other.lup);
						break;
					case MatrixType.zero, MatrixType.band1, MatrixType.upperTri, MatrixType.lowerTri:
						break;
					}
				}
				else{
					attrInitialize();
				}
			}

			/************************
			 * Initialize by a number of type "T"
			 * Params:
			 * 	num= a number
			 ************************/
			this(in T num) @nogc{
				v[]= num;
				static if(isSquare && Size > 1){
					attr.isInvertible= Ternary.no;
				}
				else{
					attr.isInvertible= isSquare? Ternary.unknown: Ternary.no;
				}
				attr.isUnitary= Ternary.no;
			}

			/************************
			 * Initialize by a static array
			 * Params:
			 * 	nums= 
			 * 	attrGiven= 
			 ************************/
			this(in TypeOfInternalArray nums, in typeof(attr) attrGiven) @nogc{
				nums[].copy(this.v[]);
				attr= attrGiven;
			}

			/// ditto
			this(in TypeOfInternalArray nums) @nogc{
			  nums[].copy(this.v[]);
			  attrInitialize();
			}

			/************************
			 * Initialize by a dynamic array
			 * Params:
			 *
			 * Throws:
			 * 	AssertError
			 ************************/
			this(in T[] nums)
		  in(nums.length == v.length, "Mismatch array length."){
			  nums[].copy(v[]);
			  attrInitialize();
			}

			/************************
			 * for only dense matrix
			 ************************/
			static if(Shape is MatrixType.dense){
				static if(MatOdr is MajorOrder.row){
					this(in T[Column][Row] num) @nogc{
					  size_t st, en;
						foreach(i; 0u..Row){
							st= i*Column;
							en= st+Column;
							num[i][].copy(this.v[st..en]);
						}
						attrInitialize();
					}

					this(in Vector!Row[Column] vecArray) @nogc{
						size_t idx= void;
						foreach(idx_h; 0..Column){
							idx= idx_h;
							foreach(idx_v; 0u..Row){
								this.v[idx]= vecArray[idx_h]._values[idx_v];
								idx += Column;
							}
						}
						attrInitialize();
					}
				}
				else static if(MatOdr is MajorOrder.column){
					this(in T[Row][Column] num) @nogc{
						size_t st, en;
						foreach(i; 0u..Column){
							st= i*Row;
							en= st+Row;
							num[i][].copy(this.v[st..en]);
						}
						attrInitialize();
					}

					this(in Vector!Row[Column] vecArray) @nogc{
						size_t idx_st= 0u, idx_en= Row;
						foreach(idx_h; 0u..Column){
							vecArray[idx_h]._values[].copy(this.v[idx_st..idx_en]);
							idx_st= idx_en;
							idx_en += Row;
						}
						attrInitialize();
					}
				}
			}	// end of block "Shape is MatrixType.dense"
		}

		/******************************************
		 * Operators & Reserved methods
		 ******************************************/
		/// plus or minus sign
		@safe pure nothrow const{
			TypeOfThis opUnary(string Op: "+")() @nogc{
				return typeof(return)(this.v[]);
			}

			TypeOfThis opUnary(string Op: "-")() @nogc{
				TypeOfInternalArray temp= void;
				temp[]= -this.v[];
				return typeof(return)(temp);
			}

			// add & sub (the same shape)
		  TypeOfThis opBinary(string Op)(in TypeOfThis rhs) @nogc
			if(isPlusOrMinusSign!Op){
			  TypeOfInternalArray temp= this.v[]; 
				mixin("temp[] " ~Op ~"= rhs.v[];");
				return typeof(return)(temp);
			}

			// add & sub (different shape)
		  TypeOfThis opBinary(string Op,
													MatrixType ShapeR: MatrixType.zero
													)(in Matrix!(Row, Column, ShapeR, MatOdr) rhs) @nogc
			if(isPlusOrMinusSign!Op && ShapeR != Shape){
				return typeof(return)();
			}

			// non-zero matrix pm zero matrix
		  Matrix!(Row, Column,
							MatOpReturnType!(Op, Shape, ShapeR),
							MatOdr) opBinary(string Op,
															 MatrixType ShapeR)(in Matrix!(Row, Column,
																														 ShapeR,
																														 MatOdr) rhs) @nogc
			if(isPlusOrMinusSign!Op && ShapeR != Shape){
				pragma(msg, "add sub from matrix.d");
				return MethodsImpl.opBinary!(PM, ShapeR)(rhs);
			}

			/// multiplication & division operators
			// matrix * scalar or matrix / scalar
			// O(m*n)
		  TypeOfThis opBinary(string Op)(in T rhs) @nogc
			if(Op == "*" || Op == "/"){
			  TypeOfInternalArray num= this.v[]; 
				mixin("num[] " ~Op ~"= rhs;");
				return typeof(return)(num);
			}

			// scalar * matrix
			// O(m*n)
		  TypeOfThis opBinaryRight(string Op: "*")(in T lhs) @nogc{
				return this.opBinary!Op(lhs);
			}

			// matrix * vector
			Vector!Row opBinary(string OP: "*")(in Vector!Column rhs) @nogc{
				auto colmat= Matrix!(Column, 1u, MatrixType.dense, MatOdr)(rhs.a);
			  return typeof(return)(this*colmat);
			}

			// matrix multiplication
			Matrix!(Row, ColumnR,
							MatrixType.zero,
							MatOdr) opBinary(string Op: "*",
															 TypeOfSize ColumnR,
															 MatrixType ShapeR: MatrixType.zero
															 )(in Matrix!(Column, ColumnR,
																						ShapeR,
																						MatOdr) rhs) @nogc{
				return typeof(return)();
			}

			Matrix!(Row, ColumnR,
							MatOpReturnType!(OP, Shape, ShapeR),
							MatOdr) opBinary(string Op: "*",
															 TypeOfSize ColumnR,
															 MatrixType ShapeR)(in Matrix!(Column, ColumnR,
																														 ShapeR, MatOdr) rhs) @nogc{
				return MethodsImpl.opBinary!(Op, ColumnR, ShapeR, MatOdr)(rhs);
			}

			/****************
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
			 ****************/
			Vector!Row opIndex(in size_t j) @nogc{
				import core.exception: onRangeError;
				typeof(return) vec;

				if(j > 0u && j <= Column){
					vec= typeof(return)(sliceColumn(j-1u));
				}
				else{
					onRangeError;
				}

				return vec;
			}

			/************************
			 * element
			 * Throws:
			 *  RangeError
			 ************************/
			T opIndex(in size_t idxRow, in size_t idxColumn) @nogc{
				typeof(return) result;
				const auto pos= MatrixPosition!(Row, Column)(idxRow, idxColumn);
				if(pos.rangeCheck){
				  auto map= indexOfInternalArray!(Row, Column, Shape, MatOdr)(pos);
					result= map.isZero? VALUE_ZERO : v[map.index];
				}
				return result;
			}

			/************************
			 * Cast operators
			 ************************/
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

			static if(MatOdr is MajorOrder.row){
				Array2D opCast(Array2D: T[Column][Row])() @nogc{
					typeof(return) num= void;

					foreach(i; 0u..Row)
						foreach(j; 0u..Column){{
								num[i][j]= this.opIndex(i+1u, j+1u);
			  	}}

					return num;
				}
			}
			else static if(MatOdr is MajorOrder.column){
				Array2D opCast(Array2D: T[Row][Column])() @nogc{
					typeof(return) num= void;
					foreach(j; 0u..Column){
						foreach(i; 0u..Row){
							num[j][i]= this.opIndex(i+1u, j+1u);
						}
					}
					return num;
				}
			}
		}	// end of block `@safe pure nothrow const'

		/******************************************
		 * Reserved methods
		 ******************************************/
		@safe const{
			Matrix!(Row, Column, Shape) dup() pure nothrow{
				return new typeof(return)(this.v);
			}

			void toString(Writer, Char)(scope Writer w, FormatSpec!Char formatSpec)
			if(isOutputRange!(Writer, const(Char)[])){
				import std.format: formatValue;
				import std.math: signbit;
				import std.array: array;
				import std.range: put, repeat;
				import std.typecons: Tuple;

			  immutable string LEFT_PAREN= "Matrix[";
			  immutable string TAB= array('\x20'.repeat(LEFT_PAREN.length));
				Tuple!(TypeOfIndex, "index", bool, "isZero") map;
				w.put(LEFT_PAREN);
				foreach(TypeOfSize i; 1u..Row+1u){
					if(i > 1u) put(w, TAB);
					w.put("[");
					foreach(TypeOfSize j; 1u..Column+1u){
						map= indexOfInternalArray!TmplArgsOfThis(MatrixPosition!(Row, Column)(i, j));
						if(map.isZero){
							w.put("*");	// always zero
						}
					  else{
							w.formatValue(v[map.index], formatSpec);
						}
						if(j < Column) w.put(", ");
					}
					w.put("]");
					if(i < Row) w.put("\n");
				}
				w.put("]");
			}

			override string toString(){
				char[] buf;
				buf.reserve(6u*Row*Column+(8u+2u*Row)*char.sizeof);
				auto fmt = FormatSpec!char("%s");
				this.toString((const(char)[] s){buf ~= s;}, fmt);
				return trustedAssumeUnique(buf);
			}
		}

		/******************************************
		 * Other common methods
		 ******************************************/
		@safe pure{
			/************************
			 * transpose
			 ************************/
			Matrix!(Column, Row,
							MatrixTypeOfTranspose!Shape,
							MatOdr) transpose() nothrow @nogc const @property{
				import std.range: iota, zip;
				enum TypeOfIndex LEN= arrayLength!(TmplArgsOfThis[0..3]);
				TypeOfInternalArray num= void;
				auto idxSet= IndexSetTranspose!TmplArgsOfThis();
				foreach(idxTr, idxSrc; zip(idxSet, iota(LEN))) num[idxTr]= this.v[idxSrc];
				return typeof(return)(num, attr);
			}

			/************************
			 * Decomposition
			 ************************/
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
						lup= LibLU.partialPivLU(this.v);
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

		/**************************
		 * Decomposed matrix
		 **************************/
		auto decomposedMatrix(DecompScheme Scheme: DecompScheme.singularValue, DecomposedMat Mat)()
		if(Mat is DecomposedMat.unitaryLeft ||
			 Mat is DecomposedMat.diagonal ||
			 Mat is DecomposedMat.unitaryRight){
			static if(Mat is DecomposedMat.unitaryLeft){
			  auto MatTmplParams= Tuple!(Row, Row, MatrixType.dense);
				T[arrayLength!MatTmplParams] num;
				auto matAttr= tuple!("isInvertible", "isUnitary")(Ternary.yes, Ternary.yes);
				return new Matrix!(MatTmplParams)(num, matAttr);
			}
			else static if(Mat is DecomposedMat.diagonal){
			  auto MatTmplParams= Tuple!(Row, Column, MatrixType.band1);
				T[arrayLength!MatTmplParams] num;
				return new Matrix!(MatTmplParams)(num);
			}
			else{
			  auto MatTmplParams= Tuple!(Column, Column, MatrixType.dense);
				T[arrayLength!MatTmplParams] num;
				auto matAttr= tuple!("isInvertible", "isUnitary")(Ternary.yes, Ternary.yes);
				return new Matrix!(MatTmplParams)(num, matAttr);
			}
		}

		auto decomposedMatrix(DecompScheme Scheme: DecompScheme.qr, DecomposedMat Mat)()
	  if(Mat is DecomposedMat.unitary ||
			 Mat is DecomposedMat.upperTri){
			static if(Mat is DecomposedMat.unitary){
				return new Matrix!(Row, Column,
													 MatrixType.dense
													 )(num, tuple!("isInvertible",
																				 "isUnitary")(Ternary.yes, Ternary.yes));
			}
			else{
				return new Matrix!(Row, Column, MatrixType.upperTri)(num);
			}
		}

		/******************************************
		 * Methods of square matrix
		 ******************************************/
		static if(isSquare){
			auto decomposedMatrix(DecompScheme Scheme: DecompScheme.lup,
														DecomposedMat Mat)() @safe pure nothrow const
		  if(Mat is DecomposedMat.lowerTri ||
				 Mat is DecomposedMat.upperTri ||
				 Mat is DecomposedMat.permutationLeft){
				static if(Mat is DecomposedMat.lowerTri){
					enum MatrixType ShapeOfReturn= MatrixType.lowerTri;
				}
				else static if(Mat is DecomposedMat.UpperTri){
					enum MatrixType ShapeOfReturn= MatrixType.upperTri;
				}
				else{
					enum MatrixType ShapeOfReturn= MatrixType.dense;
				}
				return new Matrix!(Size, Size, ShapeOfReturn, MatOdr)(lup.matrix!Mat);
			}

			/************************
			 * Trace
			 ************************/
			@property{
				T trace() @safe pure nothrow @nogc const{
					import std.range: dropOne;
					import numeric.pseudocmplx: assumeRealNum;
					typeof(return) result= void;

					if(ev is null){
						result= v[0u];
						static if(Size > 1u){
							auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();
							foreach(TypeOfIndex idx; idxSet.dropOne) result += v[idx];
						}
					}
					else{
						static if(isComplex!T) result= ev.trace;
						else result= assumeRealNum!(CT, Threshold)(ev.trace);
					}
					return result;
				}

				/************************
				 * Determinant
				 ************************/
				static if(Size < 4){
					T det() @safe pure nothrow @nogc const{
						if(attr.isInvertible is Ternary.no) return VALUE_ZERO;
						else return detImpl;
					}
				}
				else{
					T det() @safe pure nothrow{
						import numeric.pseudocmplx: assumeRealNum;
						typeof(return) result= void;

						if(attr.isInvertible is Ternary.no) result= VALUE_ZERO;
						else{
							if(ev is null) result= detImpl;	// if DIP 1008 is enable, @nogc 
							else{
								static if(isComplex!T) result= ev.det;
								else result= assumeRealNum!(CT, Threshold)(ev.det);
							}
						}

						return result;
					}
				}

				/************************
				 * Inverse matrix
				 ************************/
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

			/++++++++++++++++++++++
			 + Eigen value
			 +++++++++++++++++++++/
			void findEigenValues() @safe @property{//@safe pure nothrow @property{
				import std.stdio;	// DEBUG:
				if(ev is null){
					final switch(Shape){
					case MatrixType.zero:
						ev= new typeof(ev);
						break;
					case MatrixType.dense:
						static if(Size == 1u){
							ev= new typeof(ev)(v);
						}
						else static if(Size == 2u){
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
							foreach(idx; 0u..Size) temp[idx]= complex!RT(v[idx]);
							ev= new typeof(ev)(temp);
						}
						else ev= new typeof(ev)(v);
						break;
					case MatrixType.upperTri, MatrixType.lowerTri:
						CT[Size] num= void;
						auto idxSet= IndexSetDiag!TmplArgsOfThis();
						static if(is(T == RT)){
							foreach(TypeOfIndex idx; idxSet) num[idx]= complex!RT(v[idx]);
						}
						else{
							foreach(TypeOfIndex idx; idxSet) num[idx]= v[idx];
						}
						ev= new typeof(ev)(num);
					}
				}
				else{}	// NOP: The eigen values have been already obtained
			}

			/++
			 + eigen values
			 +/
			CT[Size] eigenValues() @safe {//@safe pure nothrow @property{
				if(ev is null) findEigenValues;
				return ev.values;
			}
		}	// end of block "static if(isSquare)"

	package:
		TypeOfInternalArray v;
		Tuple!(Ternary, "isInvertible", Ternary, "isUnitary") attr;
		static if(isSquare){
			enum LUdAlgorithm AlgoLU= LUdAlgorithm.doolittle;
			Eigen!(CT, Threshold, Size) ev;
		}

	public:
		static if(isSquare) LibLU.ArrayLUP lup;

  private:
		void attrInitialize() @safe pure nothrow @nogc{
			static if(isSquare){
				attr.isInvertible= Ternary.unknown;
				attr.isUnitary= Ternary.unknown;
			}
			else{
				attr.isInvertible= Ternary.no;
				attr.isUnitary= Ternary.unknown;
			}
		}

	  auto trustedAssumeUnique(U)(U t) @trusted pure nothrow @nogc const{
			import std.exception: assumeUnique;
			return assumeUnique(t);
		}

		// returns vector
		T[Row] sliceColumn(in size_t col) @safe pure nothrow const
		in(col >= 0u && col < Column){
			typeof(return) temp= void;

			final switch(Shape){
			case MatrixType.zero:
				temp[]= VALUE_ZERO;
				break;
			case MatrixType.dense:
				final switch(MatOdr){
				case MajorOrder.column:
					uint st= Row*col;
					temp[]= v[st .. st+Row];
					break;
				case MajorOrder.row, MajorOrder.diag:
					foreach(i; 0u..Row) temp[i]= this.opIndex(i+1u, col+1u);
				}
				break;
			case MatrixType.band1:
				foreach(i; 0u..Row){
					if(i == col) temp[i]= v[i];
					else temp[i]= VALUE_ZERO;
				}
				break;
			case MatrixType.band3, MatrixType.upperTri, MatrixType.lowerTri:
				foreach(i; 0u..Row) temp[i]= this.opIndex(i+1u, col+1u);
			}
			return temp;
		}
	} // end of class "Matrix!(size, size, shape)"

	/********************************************
	 * Global functions
	 ********************************************/
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

		return typeof(return)(num);
	}
}
