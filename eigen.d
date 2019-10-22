/*****************************************************************************
 * Eigen value & vector for Matrix type
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 *****************************************************************************/
module sekitk.eigen;

import std.complex;
import std.traits: isFloatingPoint, TemplateArgsOf;
import sekitk.base: TypeOfSize;
import numeric.pseudocmplx: isComplex;

/*************************************************************
 * Eigen values & vectors
 *************************************************************/
class Eigen(CT, real Threshold, TypeOfSize Size)
if(Size > 0 && isComplex!CT){
	import numeric.pseudocmplx: BaseRealType;
	import numeric.approx: approxEqualZero;

	alias RT= BaseRealType!CT;
	enum CT CT_ZERO= CT(0.0L);
	enum RT MAX_DIFF_ABS= cast(RT)Threshold;

	/**************
	 * Constructors
	 *************/
	@safe pure nothrow{
		// Default constructor
		this() @nogc{
			invertible= false;
			ev[]= CT_ZERO;
		}

		// Copy constructor
		this(in typeof(this) other) @nogc{
			this.invertible= other.invertible;
			if(invertible) ev[]= other.ev[];
		}

		// Eigen values (argument is Complex!RT[SIZE])
		this(in CT[Size] cnum) @nogc{
			//import std.algorithm.searching: any;
			ev[]= cnum[];
			invertible= true;
			foreach(elm; ev){	// FIXME: parallel
				if(approxEqualZero!(CT, MAX_DIFF_ABS)(elm)){
					invertible= false;
					break;
				}
			}
		}

		// Eigen values (argument is RT[SIZE])
		this(in RT[Size] rnum) @nogc{
			invertible= true;
			foreach(ubyte idx, RT elm; rnum){
				ev[idx]= complex!RT(elm);
				if(invertible && approxEqualZero!(RT, MAX_DIFF_ABS)(elm)) invertible= false;
			}
		}
	}	// @safe pure nothrow

	@property @safe pure const{
		// invertible or singular
		bool isInvertible() nothrow @nogc{return invertible;}

		// eigen values
		CT[Size] values() nothrow @nogc{return ev;}

		//(Vector!size)[size] vectors(){}

		// trace
		CT trace() nothrow @nogc{return traceOrDet!"+";}

		// determinant
		CT det() nothrow @nogc{
			typeof(return) result= void;
			if(invertible) result= traceOrDet!"*";
			else result= CT_ZERO;
			return result;
		}

		// jordanNormalForm
		// diagonalize
	}

private:
	CT traceOrDet(string OP)() @safe pure nothrow @nogc const{
		typeof(return) result= ev[0];
		static if(Size > 1u){
			mixin("foreach(num; ev[1u..$]) result " ~OP ~"= num;");
		}
		return result;
	}

	bool invertible;
	CT[Size] ev;
}
