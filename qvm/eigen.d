/*****************************************************************************
 * Eigen value & vector for Matrix type
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 *****************************************************************************/
module sekitk.qvm.eigen;

import std.complex;
import std.traits: isDynamicArray, isFloatingPoint, isStaticArray, TemplateArgsOf;
import sekitk.qvm.common: TypeOfSize;
import sekitk.complex.pseudo: isComplex;

/*************************************************************
 * Eigen values & vectors
 *************************************************************/
struct Eigen(CT, real Threshold, TypeOfSize Size)
if(Size > 0 && isComplex!CT){
	import sekitk.complex.pseudo: BaseRealType;
	import sekitk.approx: approxEqualZero;

	alias RT= BaseRealType!CT;
	enum CT CT_ZERO= CT(0.0L);
	enum RT MAX_DIFF_ABS= cast(RT)Threshold;

  //Constructors
	@safe pure nothrow @nogc{
		/// Copy constructor
		this(ref return const inout typeof(this) other){}

		/// Eigen values (argument is Complex!RT[SIZE])
		this(ArrayType)(in ArrayType cnum) @nogc
		if((isDynamicArray!ArrayType && is(ArrayType: CT[])) ||
			 (isStaticArray!ArrayType && is(ArrayType: CT[Size])))
		in(cnum.length == Size){
			import std.algorithm: any;
			_evalues[]= cnum[];
			_invertible= !_evalues[].any(a => approxEqualZero(a));
/+
			foreach(const elm; _evalues){
				if(approxEqualZero!(CT, MAX_DIFF_ABS)(elm)){
					_invertible= false;
					break;
				}
				else continue;
			}+/
		}
	}

	@property @safe pure const{
		// invertible or singular
		bool isInvertible() nothrow @nogc{return _invertible;}

		// eigen values
		CT[Size] values() nothrow @nogc{return _evalues;}

		//(Vector!size)[size] vectors(){}

		// trace
		CT trace() nothrow @nogc{return traceOrDet!"+";}

		// determinant
		CT det() nothrow @nogc{
			typeof(return) result= void;
			if(_invertible) result= traceOrDet!"*";
			else result= CT_ZERO;
			return result;
		}

		// jordanNormalForm
		// diagonalize
	}

private:
	CT traceOrDet(string Op)() @safe pure nothrow @nogc const{
		typeof(return) result= _evalues[0];
		static if(Size > 1u){
			mixin("foreach(num; _evalues[1u..$]) result " ~Op ~"= num;");
		}
		return result;
	}

	bool _invertible;
	CT[Size] _evalues;
}
