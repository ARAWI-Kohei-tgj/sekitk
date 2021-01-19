/******************************************************************************
 * Operations for integer types
 *
 * Authors: 新井 浩平 (Kohei ARAI), arawi_kohei_takasaki@yahoo.co.jp
 * Version: 1.0
 ******************************************************************************/
module sekitk.integers.intop;

import std.traits: isIntegral;

/*************************************************************
 * Combination
 *************************************************************/
Z combination(Z)(in T n, in T i) @safe pure nothrow @nogc
if(isIntegral!Z)
in(n >= i){
	import sekitk.integers.progression: factorial;
	return factorial!Z(n)/(factorial!T(cast(Z)(n-i))*factorial!Z(i));
}
