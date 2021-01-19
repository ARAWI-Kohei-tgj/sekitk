/**************************************************************
 * approximately equal
 *
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 2.1
 *************************************************************/
module sekitk.approx;

import std.complex;
import std.traits: isFloatingPoint, isImplicitlyConvertible, TemplateArgsOf;
import sekitk.complex.pseudo: isComplex, BaseRealType;

@safe pure nothrow @nogc{
	/*****************************
	 * approximately equal to zero
	 *
	 * |x| < epsilon_abs
	 * |x+yi| < epsilon_abs
	 *****************************/
  bool approxEqualZero(CT, real MaxDiffAbs)(in CT num)
	if(isComplex!CT || isFloatingPoint!CT){
		typeof(return) result= void;

		static if(isFloatingPoint!CT){
			import std.math: fabs;
			result= (num.fabs < MaxDiffAbs);
			}
		else{
			import std.complex: sqAbs;
			result= (num.sqAbs < MaxDiffAbs^^2);
		}

		return result;
	}
	unittest{
		assert(approxEqualZero!(Complex!float, 1.0e-3f)(complex!float(0.0)));
		assert(approxEqualZero!(double, 1.0e-3)(1.0e-6));
	}

	/*****************************
	 * approximately equal to one
	 *****************************/
	bool approxEqualOne(T, T threshold)(in T num, in T errRel)
	if(isFloatingPoint!T){
		import std.math: approxEqual;
		return approxEqual(num, 1.0L, errRel, threshold);
	}
}

