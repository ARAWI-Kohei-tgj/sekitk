module sekitk.qvm.impl.permutation;
/+
import std.traits: isFloatingPoint;
import numeric.pseudocmplx: isComplex;
import sekitk.base: TypeOfSize, MatrixType, MajorOrder, matrixConstraint;

mixin template MatrixImplPerm(T, T Threshold)
if(isFloatingPoint!T || isComplex!T){
	import sekitk.base;
+/
enum string MATRIX_IMPL_PERM= q{
	// permutation matrix
	class Matrix(TypeOfSize Row, TypeOfSize Column,
							 MatrixType Shape: MatrixType.permutation,
							 MajorOrder MatOdr= MajorOrder.row)
	if(matrixConstraint!(Row, Column, Shape)){
		import std.array: array, staticArray;
		import std.format: FormatSpec;
		import std.range: iota, isOutputRange;
		import std.traits: isDynamicArray, isStaticArray;
		import std.typecons: Ternary;
		import sekitk.base: Identity;
		import sekitk.qvm.attribute;

		enum T ZERO_VALUE= Identity!(T, "+");
		enum TypeOfSize Size= Row;
		alias TypeOfInternalArray= TypeOfSize[Size];
		alias TypeOfAttr= MatrixAttribute!(Row, Column, MatrixType.permutation);
		alias TypeOfThis= typeof(this);

		// constructors
		@safe pure nothrow @nogc{
			/// Default constructor
			this(){
				this._idx= iota(Size).staticArray!Size;
				this._isOdd= false;
			}

			///
			this(ArrayType)(scope inout ArrayType idx, scope inout Ternary isOdd= Ternary.unknown)
		  if((isDynamicArray!ArrayType && is(TypeOfSize[]: ArrayType)) ||
				 (isStaticArray!ArrayType && is(TypeOfSize[Size]: ArrayType)))
			in(idx.length == Size, "mismatch array length"){
				static if(is(ArrayType == TypeOfSize[]) ||
									is(ArrayType == TypeOfSize[Size])) this._idx[]= idx[];
				else{
					iota(Size).each!(i => this._idx[i]= cast(TypeOfSize)(idx[i]));
				}

				if(isOdd !is Ternary.unknown){
					this._isOdd= (isOdd is Ternary.yes);
				}
				else{
					this._isOdd= (in TypeOfInternalArray seq) @safe pure nothrow @nogc const{
						TypeOfSize counter= 0u;
						foreach(i; 0u..Size){
							foreach(j; i+1u..Size){
								if(seq[i] > seq[j]) ++counter;
							}
						}
						return counter%2u != 0;
					}(this._idx);
				}
			}
		}	// constructors

		// operators
		@safe pure nothrow const{
			// permutation * zero
			Matrix!(Row, ColumnR, MatrixType.zero, MatOdr)
				opBinary(string Op: "*",
								 TypeOfSize ColumnR,
								 MatrixType ShapeR: MatrixType.zero
								 )(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs){
				return new typeof(return)();
			}

			// permutation * permutation
			Matrix!(Row, ColumnR, MatrixType.permutation, MatOdr)
				opBinary(string Op: "*",
								 TypeOfSize ColumnR: Row,
								 MatrixType ShapeR: MatrixType.permutation
								 )(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs){
				TypeOfInternalArray idx= void;
				foreach(i; 0u..Size) idx[rhs._idx[i]]= this._idx[i];
				return new typeof(return)(idx);
			}

/+
 1  0  0  0  0
 0  1  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  1

 a  b
 c  d
 e  f
 g  h
 i  j
+/

			// matrix multiplications
			Matrix!(Row, ColumnR, MatrixType.dense, MatOdr)
				opBinary(string Op: "*",
								 TypeOfSize ColumnR,
								 MatrixType ShapeR)(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs)
		  if(ShapeR !is MatrixType.zero && ShapeR !is MatrixType.permutation){
				import std.traits: TemplateArgsOf;
			  enum MatrixType ShapeDest= TemplateArgsOf!(typeof(return))[2];
				T[arrayLength!(Row, ColumnR, MatrixType.dense)] num= void;

				return new typeof(return)(num);
			}


			/// cast to dense matrix
			Mat opCast(Mat: Matrix!(Row, Column,
															MatrixType.dense, MatOdr))(){
				T[Size][Size] nums= (){
					T[arrayLength!(TemplateArgsOf!Mat[0..3])] temp= ZERO_VALUE;
					nums= temp.chunk!Size;
				}();

				foreach(TypeOfSize i, idx; this._values) nums[i][idx]= T(1.0L);

				return new Mat(nums, _attr);
			}
		}

		// reserved methods
		@safe const{
			/****************************************
			 * Convert to string
			 ****************************************/
			override string toString(){
				import std.exception: assumeUnique;

				char[] buf;
				{
					enum size_t RESERVE_SIZE= (LEFT_PAREN.length	// "Matrix["
																		 +"[],\n".length	// "[],\n" or "[]]\n"
																		 +DELIM.length*(Column-1)	// ", "
																		 +char.sizeof*Column	// '0' or '1'
																		 )*Row;
					buf.reserve(RESERVE_SIZE);
				}
				auto fmt= FormatSpec!char("%s");
				this.toString((const(char)[] s){buf ~= s;}, fmt);

				return (char[] bufMutable) @trusted{return assumeUnique(bufMutable);}(buf);
			}

			/// ditto
			void toString(Writer, Char)(scope Writer w, FormatSpec!Char formatSpec)
			if(isOutputRange!(Writer, const(Char)[])){
				import std.range: put, repeat;
				immutable string TAB= array('\x20'.repeat(LEFT_PAREN.length));
				w.put(LEFT_PAREN);
				foreach(TypeOfSize i; 0u..Row){
					if(i > 0) w.put(TAB);
					w.put("[");
					foreach(TypeOfSize j; 0u..Column){
						w.put((_idx[i] == j)? "1": "0");
						if(j < Column-1u) w.put(", ");
					}
					w.put((i < Row-1u)? "],\n": "]");
				}
				w.put("]");
			}
		}

		// Other methods
		@safe pure nothrow const{
			/// transpose
			TypeOfThis transpose(){
				TypeOfSize[Size] idxTr= void;
				foreach(TypeOfSize i; 0u..Size) idxTr[_idx[i]]= i;

				return new typeof(return)(idxTr, Ternary(_isOdd));
			}

			// determinant
			byte det() @nogc{
				return _isOdd? -1: 1;
			}

			/// inverse of an unitary matrix is equal to its transposed one.
			alias inverse= transpose;
		}

	private:
		bool _isOdd;
		TypeOfInternalArray _idx;
		TypeOfAttr _attr;

		enum string LEFT_PAREN= "Matrix[";
		enum string DELIM= ", ";
	}
};
