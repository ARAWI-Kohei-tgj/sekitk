/******************************************************************************
 * Base class for sekitk.vector & sekitk.quaternion
 *
 * This module is a sub-module for module SekiTK.
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version 1.0:
 ******************************************************************************/
module sekitk.vectorbase;

import std.traits: isFloatingPoint;
import numeric.pseudocmplx: isComplex;
import sekitk.base: TypeOfSize;

version= IMPL_CLASS;
/**************************************************************
 * Implementation like a base class
 **************************************************************/
mixin template VectorBase(T, real Threshold, TypeOfSize Size)
if((isFloatingPoint!T || isComplex!T) & (Size > 0u)){
	enum string VECTOR_BASE_IMPL= q{
		import std.math: sqrt, approxEqual;
		import std.meta: templateAnd;
		import std.traits: Unqual, isDynamicArray, isStaticArray, isIntegral,
			ParameterTypeTuple;
		import numeric.approx: approxEqualZero;
		import sekitk.base: isPlusOrMinusSign;
		import sekitk.exceptions: ZeroNorm;

		import std.format: FormatSpec;
		import std.range.primitives: isOutputRange;
		import std.complex: Complex;
		import numeric.pseudocmplx: BaseRealType;

		static if(isComplex!T){
			alias CT= T;
			alias RT= BaseRealType!T;
		}
		else{
			alias CT= Complex!T;
			alias RT= T;
		}

		alias TypeOfThis= typeof(this);

public:
		// Constructors
		@safe pure nothrow @nogc{
			/// default constructor
			this(){}

			/// copy constructor
			this(inout scope TypeOfThis other){
				this._values[]= other._values[];
			}

			/****************************************
			 * Params:
			 *	num= a number of type T
			 ****************************************/
			this(Scalar)(in Scalar num)
			if(is(Scalar: T)){
				this._values[]= num;
			}

			/****************************************
			 * Params:
			 *	num= an array
			 ****************************************/
			this(ArrayType)(in ArrayType nums)
			if((isDynamicArray!ArrayType && is(ArrayType: T[])) ||
				 (isStaticArray!ArrayType && is(ArrayType: T[Size])))
		  in(nums.length == Size){
					this._values[]= nums[];
		  }
		}

		// Operators & reserved methods
		@safe pure nothrow @nogc{
			///
			void opOpAssign(string Op)(in TypeOfThis rhs)
			if(isPlusOrMinusSign!Op){
				mixin("this._values[] "~ Op ~"= rhs._values[];");
			}

			///
			void opOpAssign(string Op, TypeR)(in TypeR rhs)
			if((Op == "*" || Op == "/") && 
				 is(TypeR: T)){
				mixin("this._values[] "~ Op~ "= rhs;");
			}
		}

		@safe pure nothrow const{
			///
			TypeOfThis opUnary(string Op)()
			if(isPlusOrMinusSign!Op){
				typeof(return) result;
				static if(Op =="+"){
					result= new TypeOfThis(this);
				}
				else{
					T[Size] num= _values[];
					foreach(ref elm; num) elm= -elm;
					result= new TypeOfThis(num);
				}
				return result;
			}

			/****************************************
			 * Params:
			 *	rhs= other
			 *
			 * Returns:
			 *	add or sub
			 ****************************************/
			TypeOfThis opBinary(string Op, TypeR: TypeOfThis)(in TypeR rhs)
			if(isPlusOrMinusSign!Op){
				auto result= new TypeOfThis(this);
				result.opOpAssign!Op(rhs);
				return result;
			}

			/****************************************
			 * Params:
			 *	rhs= a scalar
			 ****************************************/
			TypeOfThis opBinary(string Op, TypeR)(in TypeR rhs)
			if((Op == "*" || Op == "/") && is(TypeR: T)){
				auto result= new TypeOfThis(this._values);
				result.opOpAssign!Op(rhs);
				return result;
			}

			/****************************************
			 * Params:
			 *	lhs= a scalar
			 ****************************************/
			TypeOfThis opBinaryRight(string Op: "*", TypeL)(in TypeL lhs)
			if(is(TypeL: T)){
				auto result= new TypeOfThis(this._values);
				result.opBinary!Op(lhs);
				return result;
			}

			/****************************************
			 * equality
			 *
			 * Params:
			 *	rhs= other
			 ****************************************/
			version(IMPL_CLASS){
				override bool opEquals(Object o) @nogc{
					import std.traits: ConstOf;

					if(auto rhs= cast(ConstOf!TypeOfThis)o){
						foreach(idx; 0u..Size){
							if(!approxEqualZero!(T, Threshold)(this._values[idx]-rhs._values[idx])
								 ) return false;
							else continue;
						}
						return true;
					}
					else{
						return false;
					}
				}
			}
			else{
				bool opEquals()(auto ref const TypeOfThis rhs) @nogc{
					foreach(idx; 0u..Size){
						if(!approxEqualZero!(T, Threshold)(this._values[idx]-rhs._values[idx])) return false;
						else continue;
					}
					return true;
				}
			}
		}

		@safe{
			///
			int opApply(Dg)(scope Dg dg)
			if(ParameterTypeTuple!Dg.length == 1){
				typeof(return) result= 0;
				foreach(elm; _values) result= dg(elm);

				return result;
			}

			///
			int opApplyReverse(Dg)(scope Dg dg) if(ParameterTypeTuple!Dg.length == 1){
				typeof(return) result= 0;
				foreach_reverse(elm; _values) result= dg(elm);

				return result;
			}
		}

		// Other methods
		@safe pure const{
			/****************************************
			 * Returns:
			 *	Euclidean norm
			 ****************************************/
			RT norm() nothrow @nogc @property{
				return sqrt(this.dotProdImpl(this));
			}

			/****************************************
			 * Returns:
			 *	unit vector
			 *
			 * Throws:
			 *	ZeroNorm
			 ****************************************/
			TypeOfThis unit() @property{
				T vecLen= this.norm;
				typeof(return) result= new TypeOfThis(this);

				if(!approxEqualZero!(T, Threshold)(vecLen)){
					result /= vecLen;
				}
				else{
					throw new ZeroNorm;
				}

				return result;
			}

			/****************************************
			 * Params:
			 *	rhs= other
			 ****************************************/
			bool isApproxEqualTo(in TypeOfThis rhs) nothrow @nogc{
				typeof(return) result= true;
				foreach(i; 0u..Size){
					if( !approxEqualZero!(T, Threshold)(this._values[i]-rhs._values[i]) ){
						result= false;
						break;
					}
				}
				return result;
			}

			/****************************************
			 * Returns:
			 *	true
			 ****************************************/
			bool isNormalized() nothrow @nogc @property{
				return approxEqual(this.norm, RT(1.0L));
			}
		}

	private:
		/******************************************
		 * Dot product
		 ******************************************/
		RT dotProdImpl(in TypeOfThis rhs) @safe pure nothrow @nogc const{
			import numeric.pseudocmplx: conj;
			RT sum= RT(0.0L);
			foreach(i; 0u..Size) sum += this._values[i]*rhs._values[i].conj;
			return sum;
		}

	package:
		T[Size] _values;
	};
}
