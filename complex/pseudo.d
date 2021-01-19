module sekitk.complex.pseudo;

import std.traits: isFloatingPoint;

/**************************************************************
 * Traits for complex type defined in std.complex.d
 *
 * This code is based on a code in following site.
 * See_Also: https://wiki.dlang.org/Is_expression
 **************************************************************/
template isComplex(T){
	import std.complex;
	static if(is(T == Complex!double)) enum bool isComplex= true;
	else static if(is(T == Complex!float)) enum bool isComplex= true;
	else static if(is(T == Complex!real)) enum bool isComplex= true;
	else enum bool isComplex= false;
}

/**************************************************************
 * R
 *
 * Example:
 * ---
 * assert(is(BaseRealType!(Complex!float) == float));
 * assert(is(BaseRealType!double == double));
 * ---
 **************************************************************/
template BaseRealType(T) if(isComplex!T || isFloatingPoint!T){
	static if(isComplex!T){
		import std.traits: TemplateArgsOf;
		alias BaseRealType= TemplateArgsOf!T[0];
	}
	else{
		alias BaseRealType= T;
	}
}

/**********************************************
 * Complex!CT(num, 0.0) -> num with type RT
 **********************************************/
BaseRealType!CT assumeRealNum(CT, real Threshold)(in CT num) @safe pure nothrow @nogc
if(isComplex!CT){
	import numeric.approx: approxEqualZero;
	assert(approxEqualZero!(typeof(return), Threshold)(num.im), "assumeRealNum failed: the imaginary part is not equal to zero.");
	return num.re;
}

/**************************************************************
 * Pseudo methods of complex number for real and integer number
 **************************************************************/
@safe pure nothrow @nogc @property{
	T re(T)(in T num) if(isFloatingPoint!T ||isIntegral!T){return num;}
	T im(T)(in T num) if(isFloatingPoint!T || isIntegral!T){return T(0);}
	T conj(T)(in T num) if(isFloatingPoint!T || isIntegral!T){return num;}
	T sqAbs(T)(in T num) if(isFloatingPoint!T || isIntegral!T){return num^^2;}
}
