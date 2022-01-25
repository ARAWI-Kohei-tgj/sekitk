/******************************************************************************
 * Sub-module for numeric.sekitk.qvm
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 ******************************************************************************/
module sekitk.qvm.quaternion;

import std.traits: isFloatingPoint;

mixin template QuaternionImpl(T, T Threshold)
if(isFloatingPoint!T){
  import sekitk.qvm.common: TypeOfSize;
  import sekitk.qvm.impl.vectorbase;

  /************************************************************
   * Quaternion type
   ************************************************************/
  struct Quaternion{
    mixin(VECTOR_BASE_IMPL);

  private:
    enum TypeOfSize Size= 4u;

  public:
    // constructors
    @safe pure nothrow @nogc{
      /****************************************
       * Initialize with each elements
       *
       * Params:
       *	w= scalar part
       *	x= coefficient of base i
       *	y= coefficient of base j
       *	z= coefficient of base k
       ****************************************/
      this(in T w, in T x, in T y, in T z){
	this._cmpn[0]= w;
	this._cmpn[1]= x;
	this._cmpn[2]= y;
	this._cmpn[3]= z;
      }

      /****************************************
       * Initialize by a scalar and a vector
       *
       * Params:
       *	w= scalar part
       *	vec= vector part
       ****************************************/
      this(in T w, in Vector!3u vec){
	this._cmpn[0]= w;
	this._cmpn[1]= vec._cmpn[0];
	this._cmpn[2]= vec._cmpn[1];
	this._cmpn[3]= vec._cmpn[2];
      }
    }

    // Operators & reserved methods
    @safe pure nothrow @nogc{
      /// quaternion *= quaternion or quaternion /= quaternion
      void opOpAssign(string Op)(in TypeOfThis rhs)
      if(Op == "*" || Op = "/"){
	this= this.opBinary!Op(rhs);
      }

      /// Hamilton product
      TypeOfThis opBinary(string Op: "*")(in TypeOfThis rhs) const{
	T[4] num= void;
	num[0]= this._cmpn[0]*rhs._cmpn[0]
	  -this._cmpn[1]*rhs._cmpn[1]
	  -this._cmpn[2]*rhs._cmpn[2]
	  -this._cmpn[3]*rhs._cmpn[3];
	num[1]= this._cmpn[0]*rhs._cmpn[1]
	  +this._cmpn[1]*rhs._cmpn[0]
	  +this._cmpn[2]*rhs._cmpn[3]
	  -this._cmpn[3]*rhs._cmpn[2];
	num[2]= this._cmpn[0]*rhs._cmpn[2]
	  -this._cmpn[1]*rhs._cmpn[3]
	  +this._cmpn[2]*rhs._cmpn[0]
	  +this._cmpn[3]*rhs._cmpn[1];
	num[3]= this._cmpn[0]*rhs._cmpn[3]
	  +this._cmpn[1]*rhs._cmpn[2]
	  -this._cmpn[2]*rhs._cmpn[1]
	  +this._cmpn[3]*rhs._cmpn[0];
	return typeof(return)(num);
      }

      /// quaternion / quaternion
      TypeOfThis opBinary(string Op: "/")(in TypeOfThis rhs) const{
	T[4] num= void;
	T k= rhs.norm;
	num[0]= this.dotProdImpl(rhs);
	num[1]= rhs._cmpn[0]*this._cmpn[1]
	  -rhs._cmpn[1]*this._cmpn[0]
	  -rhs._cmpn[2]*this._cmpn[3]
	  +rhs._cmpn[3]*this._cmpn[2];
	num[2]= rhs._cmpn[0]*this._cmpn[2]
	  +rhs._cmpn[1]*this._cmpn[3]
	  -rhs._cmpn[2]*this._cmpn[0]
	  -rhs._cmpn[3]*this._cmpn[1];
	num[3]= rhs._cmpn[0]*this._cmpn[3]
	  -rhs._cmpn[1]*this._cmpn[2]
	  +rhs._cmpn[2]*this._cmpn[1]
	  -rhs._cmpn[3]*this._cmpn[0];
	foreach(ref elm; num) elm /= k;
	return typeof(return)(num);
      }
    }

    // Reserved methods
    @safe const{
      /****************************************
       * convert to string
       ****************************************/
      void toString(Writer, Char)(scope Writer wrt, scope const ref FormatSpec!Char formatSpec)
      if(isOutputRange!(Writer, const(Char)[])){
	import std.format: formatValue;
	import std.range.primitives: put;

	wrt.put(PAREN_START);
	wrt.formatValue(_cmpn[0], formatSpec);
	wrt.put("; ");
	wrt.formatValue(_cmpn[1], formatSpec);
	wrt.put("i, ");
	wrt.formatValue(_cmpn[2], formatSpec);
	wrt.put("j, ");
	wrt.formatValue(_cmpn[3], formatSpec);
	wrt.put("k" ~PAREN_END);
      }
      unittest{
	double[4] foo= [1.0, 2.0, 3.0, 4.0];
	auto bar= SekiTK!double.Quaternion(foo);
	assert(bar.toString == "Quaternion(1; 2i, 3j, 4k)");
      }

      /// ditto
      string toString(){
	import std.exception: assumeUnique;

	char[] buf;
	{
	  enum size_t RESERVE_SIZE= PAREN_START.length
	    +"; i, j, k".length
	    +"-123.45678".length*4
	    +PAREN_END.length;
	  buf.reserve(RESERVE_SIZE);
	}
	auto fmt= FormatSpec!char("%s");
	toString((const(char)[] s){buf ~= s;}, fmt);
	return (char[] bufMutable) @trusted{return assumeUnique(bufMutable);}(buf);
      }
    }

    // Other methods
    @property @safe pure nothrow @nogc const{
      /****************************************
       * Each elements
       ****************************************/
      T w(){return _cmpn[0];}
      T x(){return _cmpn[1];}
      T y(){return _cmpn[2];}
      T z(){return _cmpn[3];}

      /****************************************
       * Conjugate quaternion
       ****************************************/
      TypeOfThis conj(){
	T[4] num= void;
	num[0]= this._cmpn[0];
	num[1 .. 3] = -this._cmpn[1 .. 3];
	return typeof(return)(num);
      }

      /****************************************
       * Vector part
       ****************************************/
      Vector!3u vec(){
	import std.array: staticArray;
	return typeof(return)(_cmpn[1..4].staticArray!3u);
      }
    }

  private:
    enum string PAREN_START= "Quaternion(";
    enum string PAREN_END= ")";
  }
}
