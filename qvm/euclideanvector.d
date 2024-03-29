/******************************************************************************
 * Sub-module for numeic.sekitk
 *
 * Authors: Kohei ARAI (新井 浩平), arawi_kohei_takasaki@yahoo.co.jp
 * Version: 1.1
 ******************************************************************************/
module sekitk.qvm.euclideanvector;

import std.traits: isFloatingPoint;

import sekitk.complex.pseudo: isComplex;
import sekitk.qvm.common: TypeOfSize, MatrixType, MajorOrder;

mixin template VectorImpl(T, real Threshold)
if(isFloatingPoint!T || isComplex!T){
  import sekitk.qvm.common;
  import sekitk.qvm.impl.vectorbase;
  /************************************************************
   * Template of Eucledean vector
   ************************************************************/
  struct Vector(TypeOfSize Size)
  if(Size > 0u){
    mixin(VECTOR_BASE_IMPL);

  public:
    // Constructors
    @safe pure{
      /// Initialize by a column matrix
      this(MajorOrder MatOdr)(in Matrix!(Size, 1u,
					 MatrixType.dense,
					 MatOdr) vec) @nogc nothrow{
	this._cmpn[]= vec.v[];
      }
      /+
       this(T[] args ...)
       in(args.length == Size){
       import std.algorithm: fill;
       this._cmpns[].fill(args);
       }
       +/
    }

    // Operators
    @safe pure nothrow @nogc const{
      /// dot product
      auto opBinary(string Op: "*")(in typeof(this) rhs){
	return this.dotProdImpl(rhs);
      }

      /// cross product
      static if(Size == 1u || Size == 3u || Size == 7u){
	TypeOfThis opBinary(string Op: "%")(in TypeOfThis rhs){
	  T[Size] elm= void;
	  static if(Size == 1u){
	    elm[0]= T(0.0L);
	  }
	  else static if(Size == 3u){
	    elm[0]= this._cmpn[1]*rhs._cmpn[2]
	      -this._cmpn[2]*rhs._cmpn[1];
	    elm[1]= this._cmpn[2]*rhs._cmpn[0]
	      -this._cmpn[0]*rhs._cmpn[2];
	    elm[2]= this._cmpn[0]*rhs._cmpn[1]
	      -this._cmpn[1]*rhs._cmpn[0];
	  }
	  else{
	    elm[0]= this._cmpn[1]*rhs._cmpn[2]
	      -this._cmpn[2]*rhs._cmpn[1]
	      -this._cmpn[3]*rhs._cmpn[4]
	      +this._cmpn[4]*rhs._cmpn[3]
	      -this._cmpn[5]*rhs._cmpn[6]
	      +this._cmpn[6]*rhs._cmpn[5];
	    elm[1]= -this._cmpn[0]*rhs._cmpn[2]
	      +this._cmpn[2]*rhs._cmpn[0]
	      -this._cmpn[3]*rhs._cmpn[5]
	      +this._cmpn[4]*rhs._cmpn[6]
	      +this._cmpn[5]*rhs._cmpn[3]
	      -this._cmpn[6]*rhs._cmpn[4];
	    elm[2]= this._cmpn[0]*rhs._cmpn[1]
	      -this._cmpn[1]*rhs._cmpn[0]
	      -this._cmpn[3]*rhs._cmpn[6]
	      -this._cmpn[4]*rhs._cmpn[5]
	      +this._cmpn[5]*rhs._cmpn[4]
	      +this._cmpn[6]*rhs._cmpn[3];
	    elm[3]= this._cmpn[0]*rhs._cmpn[4]
	      +this._cmpn[1]*rhs._cmpn[5]
	      +this._cmpn[2]*rhs._cmpn[6]
	      -this._cmpn[4]*rhs._cmpn[0]
	      -this._cmpn[5]*rhs._cmpn[1]
	      -this._cmpn[6]*rhs._cmpn[2];
	    elm[4]= -this._cmpn[0]*rhs._cmpn[3]
	      -this._cmpn[1]*rhs._cmpn[6]
	      +this._cmpn[2]*rhs._cmpn[5]
	      +this._cmpn[3]*rhs._cmpn[0]
	      -this._cmpn[5]*rhs._cmpn[2]
	      +this._cmpn[6]*rhs._cmpn[1];
	    elm[5]= this._cmpn[0]*rhs._cmpn[6]
	      -this._cmpn[1]*rhs._cmpn[3]
	      -this._cmpn[2]*rhs._cmpn[4]
	      +this._cmpn[3]*rhs._cmpn[1]
	      +this._cmpn[4]*rhs._cmpn[2]
	      -this._cmpn[6]*rhs._cmpn[0];
	    elm[6]= -this._cmpn[0]*rhs._cmpn[5]
	      +this._cmpn[1]*rhs._cmpn[4]
	      -this._cmpn[2]*rhs._cmpn[3]
	      +this._cmpn[3]*rhs._cmpn[2]
	      -this._cmpn[4]*rhs._cmpn[1]
	      +this._cmpn[5]*rhs._cmpn[0];
	  }
	  return typeof(return)(elm);
	}
      }
      else{
	@disable auto opBinary(string OP: "^")(in Vector!Size rhs){}
      }

      // Comparison operator
      /****************************************
       * Index operator
       *
       * Params:
       * 	idx= index, (1 ≤ idx ≤ SIZE)
       *
       * Returns:
       * 	element of "T[SIZE] _cmpn" at (idx-1).
       *
       * Throws:
       *  RangeError
       ****************************************/
      T opIndex(IdxType)(in IdxType idx)
      if(isIntegral!IdxType)
      in(idx > 0u)
      in(idx <= Size){
	return _cmpn[idx-1u];
      }

      /****************************************
       * Slice operator
       *
       ****************************************/
      T[] opSlice(IdxType)(in IdxType st, in IdxType en)
      if(isIntegral!IdxType)
      in(st > 0u)
      in(en <= Size)
      in(en > st){
	return this._cmpn[st-1u..en-1u].dup;
      }

      /// length of the array
      alias opDollar= Size;

      /****************************************
       * Cast operator
       *
       * Returns:
       *  a column matrix
       ****************************************/
      Typ opCast(Typ: Matrix!(Size, 1u, MatrixType.dense, MatOdr), MatOdr)(){
	return Typ(this._cmpn);
      }

      Typ opCast(Typ: T[Size])(){
	import std.array: staticArray;
	return this._cmpn[].staticArray!Size;
      }
    }

    // Reserved methods
    @safe const{
      /// convert to string
      void toString(Writer, Char)(scope Writer wrt, scope const ref FormatSpec!Char formatSpec)
      if(isOutputRange!(Writer, const(Char)[])){
	import std.format.write: formatValue;
	import std.range.primitives: put;

	wrt.put(PAREN_START);
	foreach(scope TypeOfSize i, scope T elm; this._cmpn){
	  wrt.formatValue(elm, formatSpec);
	  if(i < Size-1) wrt.put(DELIM);
	}
	wrt.put(PAREN_END);
      }

      /// ditto
      string toString(){
	import std.exception: assumeUnique;

	char[] buf;
	{
	  enum size_t RESERVE_SIZE= PAREN_START.length	// Vector[
	    +"-123.45678".length*Size
	    +DELIM.length*(Size-1)
	    +PAREN_END.length;	// ]^T
	  buf.reserve(RESERVE_SIZE);
	}

	const auto fmt= FormatSpec!char("%s");
	toString((in const(char)[] s){buf ~= s;}, fmt);
	return (in char[] bufMutable) @trusted{return assumeUnique(bufMutable);}(buf);
      }
      unittest{
	double[4] temp= [1.2, -3.4, 5.6, -7.8];
	auto vec= SekiTK!double.Vector!4u(temp);
	assert(vec.toString() == "Vector[1.2, -3.4, 5.6, -7.8]^T");

	import std.format: format;
	assert(format("%.2f", vec) == "Vector[1.20, -3.40, 5.60, -7.80]^T");
	assert(format("%4.1f", vec) == "Vector[ 1.2, -3.4,  5.6, -7.8]^T");
      }
    }

    // Other methods
    @safe pure nothrow const{
      Matrix!(1u, Size,
	      MatrixType.dense,
	      MatOdr) rowVector(MajorOrder MatOdr= MajorOrder.row)(){
	return new typeof(return)(this._cmpn[]);
      }

      static if(Size == 2u || Size == 3u){
	T x() @nogc{return _cmpn[0];}
	T y() @nogc{return _cmpn[1];}
      }
      static if(Size == 3u){
	T z() @nogc{return _cmpn[2];}
      }

      /**************************************
       * Outor product
       *
       * a.opProdOuter(b)= [[a[0]*b[0], a[0]*b[1], a[0]*b[2]],
       *                    [a[1]*b[0], a[1]*b[1], a[1]*b[2]],
       *                    [a[2]*b[0], a[2]*b[1], a[2]*b[2]]];
       ****************************************/
      version(future){
	Matrix!(Size, Size, MatrixType.dense, MajorOrder.row) opProdOuter(in TypeOfThis rhs){
	  T[arrayLength!(Size, Size, MatrixType.dense)] num= void;

	  foreach(scope i; 0u..Size)
	    foreach(scope j; 0u..Size){{num[i*Size+j]= this._cmpn[i]*rhs._cmpn[j];}}

	  return new typeof(return)(num);
	}
      }

      /************************
       * Geometric product
       *
       * a.opProdGeometric(b)= ;
       ************************/
      //opProdGeometric(){}

      /+
       static if(Size == 3u){
       /**
       * a.opBinary!"%"(b)= a.crossProdOp * b
       **/
       Matrix!(3, 3, MatrixType.SkewHermitian) crossProdOp() nothrow{
       T[arrayLength!(3, 3, MatrixType.SkewHermitian)] num= [-a[2], a[1], -a[0]];
       return new typeof(return)(num);
       }

       /**
       * b.opBinary!"%"(a)= a.crossProdOpRight * b
       */
       Matrix!(3, 3, MatrixType.SkewHermitian) crossProdOpRight() nothrow{
       T[arrayLength!(3, 3, MatrixType.SkewHermitian)] num= [a[2], -a[1], a[0]];
       return new typeof(return)(num);
       }
       }
+/
    }

  private:
    enum string PAREN_START= "Vector[";
    enum string PAREN_END= "]^T";
    enum string DELIM= ", ";
  }
}
