/*****************************************************************************
 * basic mathematic functions
 *
 * Authors: 新井浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 *****************************************************************************/
module sekitk.base;

import std.traits: isIntegral, isUnsigned, isFloatingPoint;

/**********************************************
 * Additive and multiplicative identity
 **********************************************/
template Identity(T, string OP)
if(isFloatingPoint!T && (OP == "+" || OP == "*")){
	static if(OP == "+") enum T Identity= T(0.0L);
	else static if(OP == "*") enum T Identity= T(1.0L);
}


