/******************************************************************************
 * A simple Eucledean vector type
 *
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 ******************************************************************************/
module sekitk.simplevector;

import std.traits: isFloatingPoint, isIntegral, isSigned;
import sekitk.base: isPlusOrMinusSign, isSignedInt;
import numeric.pseudocmplx: isComplex;

struct Vector(T, size_t Size)
if((isFloatingPoint!T || isComplex!T || isSignedInt!T) && Size > 0u){
	import std.format: FormatSpec;
	import std.range.primitives: isOutputRange;

	alias Vec= typeof(this);

	static if(isComplex!T){
		import std.traits: TemplateArgsOf;
		alias BaseType= TemplateArgsOf!T[0];
	}
	else{
		import numeric.pseudocmplx: conj;
		alias BaseType= T;
	}

	/********************************************
	 * Constructor
	 ********************************************/
	@safe pure nothrow{
		this(in T[Size] elm) @nogc{
			this.num[]= elm[];
		}

		this(in T[] elm)
		in(elm.length == Size){
		  this.num[]= elm[];
		}

		this(ref return scope inout typeof(this) another) @nogc{}
	}

	/********************************************
	 * Non-const operators
	 ********************************************/
	@safe pure nothrow @nogc{
		ref Vec opOpAssign(string Op)(in typeof(this) rhs)
		if(isPlusOrMinusSign!Op){
			mixin("foreach(idx; 0u..Size) this.num[idx] " ~Op ~"= rhs.num[idx];");
			return this;
		}

		ref Vec opOpAssign(string Op)(in T rhs)
		if(Op == "*" || Op == "/"){
			mixin("this.num[] " ~Op ~"= rhs;");
			return this;
		}
	}

	/********************************************
	 * Const operators
	 ********************************************/
	@safe pure nothrow @nogc const{
		// +Vector
	  Vec opUnary(string Op: "+")(){
			return this;
		}

		// -Vector
		Vec opUnary(string Op: "-")(){
			T[Size] buf;
			foreach(idx; 0u..Size) buf[idx]= -this.num[idx];
			return Vec(buf);
		}

		// vector pm vector
		Vec opBinary(string Op)(in typeof(this) rhs)
		if(isPlusOrMinusSign!Op){
			typeof(return) result= this;
			mixin("result " ~Op ~"= rhs;");
			return result;
		}

		// vector md scalar
		Vec opBinary(string Op, ScalarType)(in ScalarType rhs)
		if((Op == "*" || Op == "/") &&
			 (isSignedInt!ScalarType ||
				isFloatingPoint!ScalarType ||
				isComplex!ScalarType)){
			typeof(return) result= this;
			mixin("result " ~Op ~"= rhs;");
			return result;
		}

		// scalar md vector
	  Vec opBinaryRignt(string Op: "*", ScalarType)(in ScalarType rhs)
		if(isSignedInt!ScalarType ||
			 isFloatingPoint!ScalarType ||
			 isComplex!ScalarType){
			return result*rhs;
		}

		// dot product
		BaseType opBinary(string Op: "*")(in typeof(this) rhs){
			typeof(return) result= BaseType(0.0L);
			foreach(idx; 0u..Size) result += num[idx]*rhs.num[idx].conj;

			return result;
		}

		// equality
		bool opEquals()(auto ref const Vec rhs){
			typeof(return) result= true;
			foreach(idx; 0u..Size){
				if(this.num[idx] != rhs.num[idx]){
					result= false;
					break;
				}
				else continue;
			}
			return result;
		}
	}

	@safe const{
		int opApply(in int delegate(in T) @safe dg){
			typeof(return) result= 0;

			foreach(elm; this.num){
				result= dg(elm);
				if(result) break;
			}
			return result;
		}

		/// to static array
		Typ opCast(Typ: T[Size])() pure nothrow @nogc{
			typeof(return) arr= this.num[];
			return arr;
		}

		/// to Euclidean vector in Sekitk
		Typ opCast(Typ: SekiTK!(T, Eps).Vector!Size, real Eps)() pure nothrow{
			return new typeof(return)(this.num);
		}
	}

	/********************************************
	 * Reserved methods
	 ********************************************/
	@safe const{
		string toString(){
			import sekitk.base: trustedAssumeUnique;

			char[] buf;
			buf.reserve(Size*5u+(Size-1u)*DELIM.length+(PREFIX.length +SUFFIX.length)*char.sizeof);
			auto fmt= FormatSpec!char("%s");
			this.toString((const(char)[] s){buf ~= s;}, fmt);
			return trustedAssumeUnique(buf);
		}

		/// ditto
		void toString(Writer, Char)(scope Writer w, scope const ref FormatSpec!Char formatSpec)
		if(isOutputRange!(Writer, const(Char)[])){
			import std.format: formatValue;
			import std.range.primitives: put;

			w.put(PREFIX);
			foreach(size_t idx, T elm; this.num){
				w.formatValue(elm, formatSpec);
				if(idx < Size-1u) w.put(DELIM);
			}
			w.put(SUFFIX);
		}
	}

	/********************************************
	 * Other methods
	 ********************************************/
	BaseType norm() @safe pure nothrow @nogc const{
		import std.math: sqrt;
		return sqrt(this*this);
	}

private:
	T[Size] num;

	enum string PREFIX= "Vector[";
	enum string SUFFIX= "]^T";
	enum string DELIM= ", ";
}
