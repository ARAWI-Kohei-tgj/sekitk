/******************************************************************************
 * Base class for sekitk.vector & sekitk.quaternion
 *
 * This module is a sub-module for module SekiTK.
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version 1.0:
 ******************************************************************************/
module sekitk.qvm.impl.vectorbase;

import std.traits: isFloatingPoint;
import sekitk.complex.pseudo: isComplex;
import sekitk.qvm.common: TypeOfSize;

/**************************************************************
 * Implementation like a base class
 **************************************************************/
enum string VECTOR_BASE_IMPL= q{
  import std.format: FormatSpec;
  import std.range: isOutputRange;
  import std.complex: Complex;
  import std.math: sqrt;
  import std.meta: templateAnd;
  import std.traits: Unqual, isDynamicArray, isStaticArray, isIntegral, ParameterTypeTuple;

  import sekitk.qvm.common: isPlusOrMinusSign;
  import sekitk.qvm.exception: ZeroNorm;

  import sekitk.complex.pseudo: BaseRealType;

private:
  alias TypeOfThis= typeof(this);

public:
  static if(isComplex!T){
    alias CT= T;
    alias RT= BaseRealType!T;
  }
  else{
    alias CT= Complex!T;
    alias RT= T;
  }

  // Constructors
  @safe pure nothrow @nogc{
    /// copy constructor
    this(ref return scope inout TypeOfThis other){}

    /****************************************
     * Params:
     *	num= a number of type T
     ****************************************/
    this(Scalar)(scope const Scalar num)
    if(is(Scalar: T)){
      this._cmpn[]= num;
    }

    /****************************************
     * Params:
     *	num= an array
     ****************************************/
    this(ArrayType)(scope const ArrayType nums)
    if((isDynamicArray!ArrayType && is(ArrayType: T[])) ||
       (isStaticArray!ArrayType && is(ArrayType: T[Size])))
    in(nums.length == Size, "mismatch array length"){
      this._cmpn[]= nums[];
    }
  }

  // Operators & reserved methods
  @safe pure nothrow @nogc{
    ///
    void opOpAssign(string Op)(in TypeOfThis rhs)
    if(isPlusOrMinusSign!Op){
      mixin("this._cmpn[] "~ Op ~"= rhs._cmpn[];");
    }

    ///
    void opOpAssign(string Op, TypeR)(in TypeR rhs)
    if((Op == "*" || Op == "/") && is(TypeR: T)){
      mixin("this._cmpn[] "~ Op~ "= rhs;");
    }

    ///
    TypeOfThis opUnary(string Op)() const
    if(isPlusOrMinusSign!Op){
      typeof(return) result;
      static if(Op =="+"){
	result= TypeOfThis(this);
      }
      else{
	T[Size] num= _cmpn[];
	foreach(ref elm; num) elm= -elm;
	result= TypeOfThis(num);
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
    TypeOfThis opBinary(string Op, TypeR: TypeOfThis)(in TypeR rhs) const
    if(isPlusOrMinusSign!Op){
      auto result= TypeOfThis(this._cmpn);
      result.opOpAssign!Op(rhs);
      return result;
    }

    /******************************************
     * Params:
     *	rhs= a scalar
     ******************************************/
    TypeOfThis opBinary(string Op, TypeR)(in TypeR rhs) const
    if((Op == "*" || Op == "/") && is(TypeR: T)){
      auto result= TypeOfThis(this._cmpn);
      //result.opOpAssign!(Op, TypeR)(rhs);
      result *= rhs;
      return result;
    }

    /******************************************
     * Params:
     *	lhs= a scalar
     ******************************************/
    TypeOfThis opBinaryRight(string Op: "*", TypeL)(in TypeL lhs) const
    if(is(TypeL: T)){
      auto result= TypeOfThis(this._cmpn);
      result.opBinary!(Op, TypeL)(lhs);
      return result;
    }

    /******************************************
     * equality
     *
     * Params:
     *	rhs= other
     ******************************************/
    bool opEquals()(auto ref const TypeOfThis rhs) const{
      import std.math: isClose;
      foreach(scope idx; 0u..Size){
	if(!isClose(this._cmpn[idx], rhs._cmpn[idx])) return false;
	else continue;
      }
      return true;
    }
  }

  @safe{
    ///
    int opApply(Dg)(scope Dg dg)
      if(ParameterTypeTuple!Dg.length == 1){
	typeof(return) result= 0;
	foreach(scope elm; _cmpn) result= dg(elm);
	return result;
      }

    ///
    int opApplyReverse(Dg)(scope Dg dg)
    if(ParameterTypeTuple!Dg.length == 1){
      typeof(return) result= 0;
      foreach_reverse(scope elm; _cmpn) result= dg(elm);
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

    /******************************************
     * Returns:
     *	unit vector
     *
     * Throws:
     *	ZeroNorm
     ******************************************/
    TypeOfThis unit() @property{
      const RT vecLen= this.norm;
      typeof(return) result= TypeOfThis(this);

      if(!approxEqualZero(vecLen)){
	result /= vecLen;
      }
      else{
	throw new ZeroNorm;
      }

      return result;
    }

    /******************************************
     * Params:
     *	rhs= other
     ******************************************/
    bool isApproxEqualTo(in TypeOfThis rhs) nothrow @nogc{
      typeof(return) result= true;
      foreach(scope i; 0u..Size){
	if( !this.approxEqualZero(this._cmpn[i]-rhs._cmpn[i]) ){
	  result= false;
	  break;
	}
      }
      return result;
    }

    /****************************************
     * Returns:
     *	true
     ******************************************/
    bool isNormalized(in RT errRel= 1e-3) nothrow @nogc @property{
      import sekitk.approx: approxEqualOne;
      return approxEqualOne!(RT, Threshold)(this.norm, errRel);
    }
  }

private:
  @safe pure nothrow @nogc const{
    alias approxEqualZero= (in T num){
      import sekitk.approx: approxEqualZero;
      return approxEqualZero!(T, Threshold)(num);
    };

    /******************************************
     * Dot product
     ******************************************/
    RT dotProdImpl(in TypeOfThis rhs){
      import sekitk.complex.pseudo: conj;
      RT sum= RT(0.0L);
      foreach(scope i; 0u..Size) sum += this._cmpn[i]*rhs._cmpn[i].conj;
      return sum;
    }
  }

package:
  T[Size] _cmpn;
};
