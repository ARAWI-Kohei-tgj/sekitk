 /*****************************************************************************
 * Solver of Quadratic equations
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 *****************************************************************************/
module mathematic.equations;

import std.math;
import std.complex;
import std.traits: isFloatingPoint;
import numeric.pseudocmplx: isComplex;

/*************************************************************
 * y= ax^2 +bx +c
 *************************************************************/
enum SolutionType: ubyte{OneR, RepeatedR, TwoR, OneC, RepeatedC, TwoC}

class QuadraticEq(T, real threshold) if(isFloatingPoint!T || isComplex!T){
	import std.traits: TemplateArgsOf; 
	import numeric.pseudocmplx: BaseRealType;

	static if(isFloatingPoint!T){
		alias CT= Complex!T;
	}
	else{
		alias CT= T;
	}
	alias BT= BaseRealType!T;

	enum BT MAX_DIFF_ABS= cast(BT)threshold;

	this(in T a, in T b, in T c) @safe pure nothrow @nogc{
		import numeric.approx: approxEqualZero;
	  coeff[2]= a;
		coeff[1]= b;
		coeff[0]= c;

		if(approxEqualZero!(T, MAX_DIFF_ABS)(coeff[2])){	// a=0
			static if(is(T == BT)) sol_typ= SolutionType.OneR;
			else sol_typ= SolutionType.OneC;
		}
		else{
			T d= discriminant;
			static if(is(T == BT)){
				if(approxEqualZero!(BT, MAX_DIFF_ABS)(d)) sol_typ= SolutionType.RepeatedR;
				else if(d > 0.0) sol_typ= SolutionType.TwoR;
				else sol_typ= SolutionType.TwoC;
			}
			else{
				if(approxEqualZero!(CT, MAX_DIFF_ABS)(d)) sol_typ= SolutionType.RepeatedC;
				else sol_typ= SolutionType.TwoC;
			}
		}
	}

  @safe pure nothrow @nogc const{
	  T value(in T x){return coeff[2]*x^^2 +coeff[1]*x +coeff[0];}

		CT[2] solve() @property{
			typeof(return) result= void;

			final switch(sol_typ){
			case SolutionType.OneR, SolutionType.OneC:
				static if(is(T == BT)) result[0]= complex(-coeff[0]/coeff[1]);
				else result[0]= -coeff[0]/coeff[1];
				result[1]= complex(BT.nan, BT.nan);
				break;
			case SolutionType.RepeatedR, SolutionType.RepeatedC:
				static if(is(T == BT)) result[0]= complex!BT(-0.5*coeff[1]/coeff[2]);
				else result[0]= -0.5*coeff[1]/coeff[2];
				result[1]= result[0];
				break;
			case SolutionType.TwoR:
				static if(is(T == BT)){
					immutable string eq0= "0.5*(-coeff[1] ", eq1= "d_sqrt)/coeff[2]";
					T d_sqrt= sqrt(discriminant);
					mixin("result[0]= complex!BT(" ~eq0 ~"-" ~eq1 ~");");
					mixin("result[1]= complex!BT(" ~eq0 ~"+" ~eq1 ~");");
				}
				else assert(0);
				break;
			case SolutionType.TwoC:
				static if(is(T == BT)){
					BT d_sqrt= sqrt(-discriminant);	// std.math.sqrt
					result[0]= complex!(BT, BT)(-coeff[1], -d_sqrt);
					result[1]= complex!(BT, BT)(-coeff[1], d_sqrt);
					foreach(ref num; result) num /= 2.0*coeff[2];
				}
				else{
					immutable string eq0= "0.5*(-coeff[1] ", eq1= "d_sqrt)/coeff[2];";
					CT d_sqrt= sqrt(discriminant);	// std.complex.sqrt
					mixin("result[0]= " ~eq0 ~"-" ~eq1);
					mixin("result[1]= " ~eq0 ~"+" ~eq1);
				}
			}
			return result;
		}
	}

	private:
	T discriminant() @safe pure nothrow @nogc @property const{
		return coeff[1]^^2 -4.0*coeff[2]*coeff[0];
	}

	T[3] coeff;
	SolutionType sol_typ;
}

//@safe pure nothrow
unittest{
 Test_R: {
		Complex!double[2] result= void;
		QuadraticEq!double eq;

		import std.math: sqrt, approxEqual;
		// P(x)= x^2-1/2x-5/2
		eq= new QuadraticEq!double(1.0, -0.5, -2.5);
	  result= eq.solve;
	  assert(approxEqual(result[0].re, 0.25-sqrt(41.0)/4.0) &&
					 approxEqual(result[1].re, 0.25+sqrt(41.0)/4.0));

		// P(x)= x^2-2\sqrt{2}x +2
		eq= new QuadraticEq!double(1.0, -2.0*sqrt(2.0), 2.0);
	  result= eq.solve;
	  assert(approxEqual(result[0].re, sqrt(2.0)) && result[0].im == 0.0 &&
					 result[1] == result[0]);

		// P(x)= 2x^2 +3x +4
		eq= new QuadraticEq!double(2.0, 3.0, 4.0);
		result= eq.solve;
		assert(approxEqual(result[0].re, -3.0/4.0) && approxEqual(result[0].im, -sqrt(23.0)/4.0));
		assert(approxEqual(result[1].re, -3.0/4.0) && approxEqual(result[1].im, sqrt(23.0)/4.0));
	}

 Test_C: {
		// P(z)= (1+i)z^2 +(-2+i)z +(-1+2i)
		Complex!double[2] results;
		Complex!double a, b, c;
		a= Complex!double(1.0, 1.0);
		b= Complex!double(-2.0, 1.0);
		c= Complex!double(-1.0, 2.0);
		auto eq= new QuadraticEq!(Complex!double)(a, b, c);
		results= eq.solve;
		assert(approxEqual(results[0].re, -0.5) && approxEqual(results[0].im, 0.5));
		assert(approxEqual(results[1].re, 1.0) && approxEqual(results[1].im, -2.0));
	}
}

/*************************************************************
 * y= ax +b
 *************************************************************/
class LinearEq(T) if(isFloatingPoint!T || isComplex!T){
	this(in T a, in T b) @safe pure nothrow @nogc{
		coeff[0]= b;
		coeff[1]= a;
	}

	@safe pure nothrow @nogc const{
		T value(in T x){return coeff[1]*x +coeff[0];}
		T solve() @property{return -coeff[0]/coeff[1];}
	}

private:
	T[2] coeff;
}
