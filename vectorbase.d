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

/**************************************************************
 * Implementation like a base class
 **************************************************************/
mixin template VectorBase(T, real Threshold, TypeOfSize Size)
if((isFloatingPoint!T || isComplex!T) & (Size > 0u)){
	import std.math: sqrt, approxEqual;
	import std.traits: Unqual;
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
	@safe pure nothrow{
		/******************************************
		 * Params:
		 *	num= a number of type T
		 ******************************************/
		this(in T num) @nogc{
			this._values[]= num;
		}

		/******************************************
		 * Params:
		 *	num= a static array
		 ******************************************/
		this(in T[Size] nums) @nogc{
			this._values[]= nums[];
		}

		/******************************************
		 * Params:
		 *	num= a dynamic array which length equals "Size"
		 *
		 * Throws:
		 *	AssertError
		 ******************************************/
		this(in T[] nums)
		in(nums.length == Size){
			import std.algorithm: copy;
			nums.copy(this._values[]);
	  }
	}

	// Operators & reserved methods
	@safe pure nothrow @nogc{
		///
		void opOpAssign(string Op)(in TypeOfThis rhs)
		if(isPlusOrMinusSign!Op){
			mixin("this.values[] "~ Op ~"= rhs.values[];");
		}

		///
		void opOpAssign(string Op)(in T rhs)
		if(Op == "*" || Op == "/"){
			mixin("this.values[] "~ Op~ "= rhs;");
		}
	}

	@safe pure nothrow @nogc const{
		///
		TypeOfThis opUnary(string Op)()
		if(isPlusOrMinusSign!Op){
			typeof(return) result;
			static if(Op =="+"){
				result= TypeOfThis(this._values);
			}
			else{
				result= TypeOfThis(-this._values[]);
			}
			return result;
		}

		/******************************************
		 * Params:
		 *	rhs= other
		 *
		 * Returns:
		 *	add or sub
		 ******************************************/
		TypeOfThis opBinary(string Op)(in TypeOfThis rhs)
		if(isPlusOrMinusSign!Op){
			auto result= TypeOfThis(this._values);
			result.opOpAssign!Op(rhs);
			return result;
		}

		/******************************************
		 * Params:
		 *	rhs= a scalar
		 ******************************************/
		TypeOfThis opBinary(string Op)(in T rhs)
		if(Op == "*" || Op == "/"){
			auto result= TypeOfThis(this._values);
			result.opOpAssign!Op(rhs);
			return result;
		}

		/******************************************
		 * Params:
		 *	lhs= a scalar
		 ******************************************/
		TypeOfThis opBinaryRight(string Op: "*")(in T lhs) nothrow{
			auto result= TypeOfThis(this._values);
			return result.opBinary!Op(lhs);
		}

		/******************************************
		 * equality
		 *
		 * Params:
		 *	rhs= other
		 ******************************************/
		bool opEquals()(auto ref const TypeOfThis rhs) nothrow @nogc{
			foreach(idx; 0u..Size){
				if(!approxEqualZero!(T, Threshold)(this._values[idx]-rhs._values[idx])) return false;
				else continue;
			}
			return true;
		}

		///
		int opApply(Dg)(scope Dg dg) if(ParameterTypeTuple!Dg.length == 1){
			typeof(return) result= 0;

			foreach(elm; a){
				result= dg(elm);
				if(result) break;
			}
			return result;
		}

		///
		int opApplyReverse(Dg)(scope Dg dg) if(ParameterTypeTuple!Dg.length == 1){
			typeof(return) result= 0;

			foreach_reverse(elm; a){
				result= dg(elm);
				if(result) break;
			}
			return result;
		}
	}

	// Other methods
	@safe pure const{
		/******************************************
		 * Returns:
		 *	Euclidean norm
		 ******************************************/
		RT norm() nothrow @nogc @property{
			return sqrt(this.dotProdImpl(this));
		}

		/******************************************
		 * Returns:
		 *	unit vector
		 *
		 * Throws:
		 *	ZeroNorm
		 ******************************************/
	  TypeOfThis unit() @property{
			T vecLen= this.norm;
			T[Size] temp= this._values[];

			if(!approxEqualZero!(T, Threshold)(vecLen)){
				temp[] /= vecLen;
			}
			else{
				throw new ZeroNorm;
			}

			return typeof(return)(temp);
		}

		/******************************************
		 * Params:
		 *	rhs= other
		 ******************************************/
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

		/******************************************
		 * Returns:
		 *	true
		 ******************************************/
		bool isNormalized() nothrow @nogc @property{
			return approxEqual(this.norm, RT(1.0L));
		}
	}

private:
	/********************************************
	 * Dot product
	 ********************************************/
	RT dotProdImpl(in TypeOfThis rhs) @safe pure nothrow @nogc const{
		import numeric.pseudocmplx: conj;
		RT sum= RT(0.0L);
		foreach(i; 0u..Size) sum += this._values[i]*rhs._values[i].conj;
		return sum;
	}

package:
	T[Size] _values;
}
