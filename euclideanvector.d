/******************************************************************************
 * Sub-module for numeic.sekitk
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.1
 ******************************************************************************/
module sekitk.euclideanvector;

import std.traits: isFloatingPoint;

import numeric.pseudocmplx: isComplex;
import sekitk.base: TypeOfSize, MatrixType, MajorOrder;

mixin template VectorImpl(T, real Threshold)
if(isFloatingPoint!T || isComplex!T){
	import sekitk.base;
	import sekitk.vectorbase: VectorBase;
	/************************************************************
	 * Template of Eucledean vector
	 ************************************************************/
	class Vector(TypeOfSize Size) if(Size > 0u){
		mixin VectorBase!(T, Threshold, Size);
		mixin(VECTOR_BASE_IMPL);

		// Constructors
		nothrow{
			import core.vararg;//: _arguments, va_arg, _argptr;

			/// Initialize by a column matrix
			this(MajorOrder MatOdr)(in Matrix!(Size, 1u,
																				 MatrixType.dense,
																				 MatOdr) vec) @safe pure @nogc{
				this._values[]= vec.v[];
			}
/+
			this(...)	@system	// FIXME:
			in(_arguments.length == Size){
				foreach(ref T num; this._values) num= va_arg!T(_argptr);
			}
+/
		}

		// Operators
	  @safe pure nothrow const{
			/// dot product
			auto opBinary(string Op: "*")(in typeof(this) rhs) @nogc{
				return this.dotProdImpl(rhs);
			}

			/// cross product
			static if(Size == 1u || Size == 3u || Size == 7u){
				TypeOfThis opBinary(string Op: "%")(in TypeOfThis rhs){
					T[Size] elm= void;
					static if(Size == 1u){
						elm[0]= T(0.0L);
					}
					else static if(Size == 3u){
						elm[0]= this._values[1]*rhs._values[2]
							-this._values[2]*rhs._values[1];
						elm[1]= this._values[2]*rhs._values[0]
							-this._values[0]*rhs._values[2];
						elm[2]= this._values[0]*rhs._values[1]
							-this._values[1]*rhs._values[0];
					}
					else{
						elm[0]= this._values[1]*rhs._values[2]
							-this._values[2]*rhs._values[1]
							-this._values[3]*rhs._values[4]
							+this._values[4]*rhs._values[3]
							-this._values[5]*rhs._values[6]
							+this._values[6]*rhs._values[5];
						elm[1]= -this._values[0]*rhs._values[2]
							+this._values[2]*rhs._values[0]
							-this._values[3]*rhs._values[5]
							+this._values[4]*rhs._values[6]
							+this._values[5]*rhs._values[3]
							-this._values[6]*rhs._values[4];
						elm[2]= this._values[0]*rhs._values[1]
							-this._values[1]*rhs._values[0]
							-this._values[3]*rhs._values[6]
							-this._values[4]*rhs._values[5]
							+this._values[5]*rhs._values[4]
							+this._values[6]*rhs._values[3];
						elm[3]= this._values[0]*rhs._values[4]
							+this._values[1]*rhs._values[5]
							+this._values[2]*rhs._values[6]
							-this._values[4]*rhs._values[0]
							-this._values[5]*rhs._values[1]
							-this._values[6]*rhs._values[2];
						elm[4]= -this._values[0]*rhs._values[3]
							-this._values[1]*rhs._values[6]
							+this._values[2]*rhs._values[5]
							+this._values[3]*rhs._values[0]
							-this._values[5]*rhs._values[2]
							+this._values[6]*rhs._values[1];
						elm[5]= this._values[0]*rhs._values[6]
							-this._values[1]*rhs._values[3]
							-this._values[2]*rhs._values[4]
							+this._values[3]*rhs._values[1]
							+this._values[4]*rhs._values[2]
							-this._values[6]*rhs._values[0];
						elm[6]= -this._values[0]*rhs._values[5]
							+this._values[1]*rhs._values[4]
							-this._values[2]*rhs._values[3]
							+this._values[3]*rhs._values[2]
							-this._values[4]*rhs._values[1]
							+this._values[5]*rhs._values[0];
					}
					return new typeof(return)(elm);
				}
			}


			// wedge product
			/*
			T[Size] opBinary(string OP: "^")(in Vector!Size rhs){}
			*/

			/************************
			 * Comparison operator
			 ************************/

			/************************
			 * Index operator
			 *
			 * Params:
			 * 	idx= index, (1 ≤ idx ≤ SIZE)
			 *
			 * Returns:
			 * 	element of "T[SIZE] _values" at (idx-1).
			 *
			 * Throws:
			 *  RangeError
			 ************************/
			T opIndex(IdxType)(in IdxType idx) @nogc
			if(isIntegral!IdxType)
			in(idx > 0u)
			in(idx <= Size){
			  return _values[idx-1u];
		  }

			/************************
			 * Slice operator
			 *
			 ************************/
			T[] opSlice(IdxType)(in IdxType st, in IdxType en)
			if(isIntegral!IdxType)
		  in(st > 0u)
			in(en <= Size)
		  in(en > st){
				return this._values[st-1u..en-1u].dup;
		  }

			/// length of the array
			alias opDollar= Size;

			/************************
			 * Cast operator
			 *
			 * Returns:
			 *  a column matrix
			 ************************/
		  Typ opCast(Typ: Matrix!(Size, 1u, MatrixType.dense, MatOdr), MatOdr)(){
				return Typ(this._values);
			}

			Typ opCast(Typ: T[Size])() @nogc{
			  Typ arr= this._values[];
				return arr;
			}
		}

	  // Reserved methods
		@safe const{
			void toString(Writer, Char)(scope Writer w, FormatSpec!Char formatSpec)
			if(isOutputRange!(Writer, const(Char)[])){
				import std.format: formatValue;
				import std.range.primitives: put;

				w.put(PREFIX);
				foreach(size_t idx, T elm; this._values){
					w.formatValue(elm, formatSpec);
					if(idx < Size-1u) w.put(DELIM);
				}
				w.put(SUFFIX);
			}

			/// ditto
		  override string toString(){
				import sekitk.base: trustedAssumeUnique;
				enum ubyte DIGITS= 6u;	// 暫定桁数

				char[] buf;
				buf.reserve(Size*6u+(Size-1u)*DELIM.length+(PREFIX.length+SUFFIX.length)*char.sizeof);
				auto fmt = FormatSpec!char("%s");
				this.toString((const(char)[] s){buf ~= s;}, fmt);
				return trustedAssumeUnique(buf);
			}

			unittest{
				double[4] temp= [1.2, -3.4, 5.6, -7.8];
				auto vec= SekiTK!double.Vector!4u(temp);
				assert(vec.toString() == "Vector[1.2, -3.4, 5.6, -7.8]^T");

				import std.format: format;
				assert(format("%.2f", vec) == "Vector[1.20, -3.40, 5.60, -7.80]^T");
				assert(format("%4.1f", vec) == "Vector[ 1.2, -3.4,  5.6, -7.8]^T");
			}
		}

	  // Other methods
		@safe pure nothrow const{
			Matrix!(1u, Size,
							MatrixType.dense,
							MatOdr) rowVector(MajorOrder MatOdr= MajorOrder.row)(){
				return new typeof(return)(this._values[]);
			}

			static if(Size == 2u || Size == 3u){
				T x() @nogc{return _values[0];}
				T y() @nogc{return _values[1];}
			}
			static if(Size == 3u){
				T z() @nogc{return _values[2];}
			}

			version(future){
			/************************
			 * Outor product
			 *
			 * a.opProdOuter(b)= [[a[0]*b[0], a[0]*b[1], a[0]*b[2]],
			 *                    [a[1]*b[0], a[1]*b[1], a[1]*b[2]],
			 *                    [a[2]*b[0], a[2]*b[1], a[2]*b[2]]];
			 ************************/
				Matrix!(Size, Size, MatrixType.dense, MajorOrder.row) opProdOuter(in TypeOfThis rhs){
					T[arrayLength!(Size, Size, MatrixType.dense)] num= void;

					foreach(i; 0u..Size)
						foreach(j; 0u..Size){{num[i*Size+j]= this._values[i]*rhs._values[j];}}

					return new typeof(return)(num);
				}
			}

			/************************
			 * Geometric product
			 *
			 * a.opProdGeometric(b)= ;
			 ************************/
			//opProdGeometric(){}

			/+
			static if(Size == 3u){
				/**
				 * a.opBinary!"%"(b)= a.crossProdOp * b
				 **/
				Matrix!(3, 3, MatrixType.SkewHermitian) crossProdOp() nothrow{
					T[arrayLength!(3, 3, MatrixType.SkewHermitian)] num= [-a[2], a[1], -a[0]];
					return new typeof(return)(num);
				}

				/**
			 	 * b.opBinary!"%"(a)= a.crossProdOpRight * b
				 */
				Matrix!(3, 3, MatrixType.SkewHermitian) crossProdOpRight() nothrow{
					T[arrayLength!(3, 3, MatrixType.SkewHermitian)] num= [a[2], -a[1], a[0]];
					return new typeof(return)(num);
				}
			}
+/
		}

	private:
		enum string PREFIX= "Vector[";
		enum string SUFFIX= "]^T";
		enum string DELIM= ", ";
	}
}
