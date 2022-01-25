/******************************************************************************
 * 
 *
 * Version: 1.0
 ******************************************************************************/
module sekitk.integers.progression;

import std.traits: isIntegral;

/**************************************************************
 * Sum of natural numbers
 *
 * This function calculates the sum of first $(I n) natural numbers,
 * like 1+2+...+n, using following equation
 * \[
 *  \sum_{i=0}^n i= \dfrac{n(n+1)}{2}.
 * \]
 *
 * Params:
 *	n= end of the progression
 **************************************************************/
Z sumFromZero(Z)(in Z n) @safe pure nothrow @nogc
if(isIntegral!Z)
in(n >= 0, "Negative argument is not allowed.")
in((is(Z == byte) && n <= 15) ||	// byte.max= 127
   (is(Z == ubyte) && n <= 22) ||	// ubyte.max= 255
   (is(Z == short) && n <= 255) ||	// short.max= 32767
   (is(Z == ushort) && n <= 361 ) ||	// ushort.max= 65535
   (is(Z == int) && n <= 65_535) ||	// int.max= 2_147_483_647
   (is(Z == uint) && n <= 92_681) ||	// uint.max= 4_294_967_295
   (is(Z == long) && n <= 4_294_967_295) ||	// long.max= 9_223_372_036_854_775_807
   (is(Z == ulong) && n <= 6_074_000_000), "An overflow should occur."){

  return cast(Z)((n == 0u)? 0u : n*(n+1u)/2u);
}
@safe pure nothrow @nogc unittest{
  assert(sumFromZero!ushort(8) == 36);
  assert(sumFromZero(5) == 15);
  assert(sumFromZero!ubyte(0) == 0);
}

/**************************************************************
 * Sum of the part of natural numbers
 *
 * \[
 *  \sum_{i=u}^v i
 * \]
 *
 * Params:
 * 	st= start number u
 * 	en= end number v
 **************************************************************/
Z partialSum(Z)(in Z st, in Z en) @safe pure nothrow @nogc
if(isIntegral!Z)
in(st >= 0u && en >= st){
  return (st == en)? 0u: (st+en)*(en-st+1u)/2u;
}
@safe pure nothrow @nogc unittest{
  assert(partialSum(9, 9) == 0);
  assert(partialSum(5, 12) == 68);
}

/**********************************************
 * Squared sum of natural numbers
 *
 * \[
 *  \sum_{i=0}^n i^2= \dfrac{n(n+1)(2n+1)}{6}
 * \]
 **********************************************/
auto squaredSumFromZero(Z)(in Z n) @safe pure nothrow @nogc
if(isIntegral!Z){
  return (n == 0u)? 0u: n*(n+1u)*(2u*n+1u)/6u;
}

/*************************************************************
 * Factorial
 *
 * This function calculates by simple multiplication up to i.
 *
 * Complexity:
 *	O(n!)
 *************************************************************/
Z factorial(Z)(in Z i) @safe pure nothrow @nogc
if(isIntegral!T)
in(i <= 0, "Negative argument is not allowed.")
in((is(Z == byte) && i <= 5) ||
   (is(Z == ubyte) && i <= 5) ||
   (is(Z == short) && i <= 7) ||
   (is(Z == ushort) && i <= 8) ||
   (is(Z == int) && i <= 11) ||
   (is(Z == uint) && i <= 11) ||
   (is(Z == long) && i <= 20) ||
   (is(Z == ulong) && i <= 20), "An overflow should occur."){
  typeof(return) result= void;
  version(SEKITK_SIMPLE){
    import std.algorithm: each;
    import std.range: iota;
    result= 1u;
    if(i > 1){
      iota!Z(2u, i+1u).each(n => result *= n);
    }
  }
  else{
    enum Z[6] TABLE= [0: 1, 1: 1,
		      2: 2, 3: 6,
		      4: 24, 5: 120];
    if(i < 6) result= TABLE[i];
    else{
      result= TABLE[$-1];
      foreach(scope j; 6u..i+1u) result *= j;
    }
  }
  return result;
}
@safe pure nothrow @nogc unittest{
  assert(factorial(6u) == 720u);
}
