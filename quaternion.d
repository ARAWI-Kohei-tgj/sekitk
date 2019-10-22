/******************************************************************************
 * Sub-module for numeric.sekitk.qvm
 *
 * Authors: 新井 浩平 (Kohei ARAI), arai-kohei-xg@ynu.jp
 * Version: 1.0
 ******************************************************************************/
module sekitk.quaternion;

import std.traits: isFloatingPoint;

mixin template QuaternionImpl(T, T Threshold) if(isFloatingPoint!T){
	import sekitk.base: TypeOfSize;
	import sekitk.vectorbase: VectorBase;

	/************************************************************
	 * Quaternion type
	 ************************************************************/
	struct Quaternion{
	private:
		enum TypeOfSize Size= 4u;

	public:
		mixin VectorBase!(T, Threshold, 4u);

		// constructors
		@safe pure nothrow @nogc{
			/****************************************
			 * Initialize by each elements
			 *
			 * Params:
			 *	w= scalar part
			 *	x= coefficient of base i
			 *	y= coefficient of base j
			 *	z= coefficient of base k
			 ****************************************/
			this(in T w, in T x, in T y, in T z){
				this._values[0]= w;
				this._values[1]= x;
				this._values[2]= y;
				this._values[3]= z;
			}

			/****************************************
			 * Initialize by a scalar and a vector
			 *
			 * Params:
			 *	w= scalar part
			 *	vec= vector part
			 ****************************************/
			this(in T w, in Vector!3u vec){
				this._values[0]= w;
				this._values[1]= vec._values[0];
				this._values[2]= vec._values[1];
				this._values[3]= vec._values[2];
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
				num[0]= this._values[0]*rhs._values[0]
					-this._values[1]*rhs._values[1]
					-this._values[2]*rhs._values[2]
					-this._values[3]*rhs._values[3];
				num[1]= this._values[0]*rhs._values[1]
					+this._values[1]*rhs._values[0]
					+this._values[2]*rhs._values[3]
					-this._values[3]*rhs._values[2];
				num[2]= this._values[0]*rhs._values[2]
					-this._values[1]*rhs._values[3]
					+this._values[2]*rhs._values[0]
					+this._values[3]*rhs._values[1];
				num[3]= this._values[0]*rhs._values[3]
					+this._values[1]*rhs._values[2]
					-this._values[2]*rhs._values[1]
					+this._values[3]*rhs._values[0];
				return typeof(return)(num);
			}

			/// quaternion / quaternion
			TypeOfThis opBinary(string Op: "/")(in TypeOfThis rhs) const{
				T[4] num= void;
				T k= rhs.norm;
				num[0]= this.dotProdImpl(rhs);
				num[1]= rhs._values[0]*this._values[1]
					-rhs._values[1]*this._values[0]
					-rhs._values[2]*this._values[3]
					+rhs._values[3]*this._values[2];
				num[2]= rhs._values[0]*this._values[2]
					+rhs._values[1]*this._values[3]
					-rhs._values[2]*this._values[0]
					-rhs._values[3]*this._values[1];
				num[3]= rhs._values[0]*this._values[3]
					-rhs._values[1]*this._values[2]
					+rhs._values[2]*this._values[1]
					-rhs._values[3]*this._values[0];
				foreach(ref elm; num) elm /= k;
				return typeof(return)(num);
			}
		}

		// Reserved methods
		@safe const{
			TypeOfThis dup() pure nothrow{
				return typeof(return)(this._values);
			}

			/****************************************
			 * 
			 ****************************************/
			void toString(Writer, Char)(scope Writer wrt, FormatSpec!Char formatSpec)
			if(isOutputRange!(Writer, const(Char)[])){
				import std.format: formatValue;
				import std.range.primitives: put;

				uint i;	// =0
				put(wrt, PREFIX);
				formatValue(wrt, _values[0], formatSpec);
				put(wrt, "; ");
				formatValue(wrt, _values[1], formatSpec);
				put(wrt, "i, ");
				formatValue(wrt, _values[2], formatSpec);
				put(wrt, "j, ");
				formatValue(wrt, _values[3], formatSpec);
				put(wrt, "k" ~SUFFIX);
			}

			/// ditto
			string toString(){
				import numeric.sekitk.qvm_base: trustedAssumeUnique;

				char[] buf;
				buf.reserve(Size*6u+(PREFIX.length+SUFFIX.length)*char.sizeof);
				auto fmt = FormatSpec!char("%s");
				this.toString((const(char)[] s){buf ~= s;}, fmt);
				return trustedAssumeUnique(buf);
			}

		  unittest{
				double[4] foo= [1.0, 2.0, 3.0, 4.0];
				auto bar= SekiTK!double.Quaternion(foo);
				assert(bar.toString == "Quaternion(1; 2i, 3j, 4k)");
			}
		}

		// Other methods
		@property @safe pure nothrow @nogc const{
			/****************************************
			 * Each elements
			 ****************************************/
			T w(){return _values[0];}
			T x(){return _values[1];}
			T y(){return _values[2];}
			T z(){return _values[3];}

			/****************************************
			 * Conjugate quaternion
			 ****************************************/
			TypeOfThis conj(){
				T[4] num= void;
				num[0]= this._values[0];
				num[1 .. 3] = -this._values[1 .. 3];
				return typeof(return)(num);
			}

			/****************************************
			 * Vector part
			 ****************************************/
			Vector!3u vec(){
				T[3] temp= _values[1..4];
				return Vector!3u(temp);
			}
		}

	private:
		enum string PREFIX= "Quaternion(";
		enum string SUFFIX= ")";
	}
}
