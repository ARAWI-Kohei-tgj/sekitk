/******************************************************************************
 * Exceptions
 ******************************************************************************/
module sekitk.qvm.exception;
import sekitk.qvm.common: TypeOfSize, MatrixType;

/**************************************************************
 * ZeroNorm
 *
 * This exception class will be thrown if property function "unit"
 * of a zero length vector is called.
 *
 * Examples:
 * ---
 * auto foo= SekiTK!double.Vector2(0.0, 0.0);
 * try{
 *     auto bar= foo.unit;
 * }
 * catch(ZeroNorm err){writeln(err.msg);}
 * ---
 **************************************************************/
class ZeroNorm: Exception{
	this() @safe pure nothrow @nogc{
		super("Zero length vector.", __FILE__, __LINE__, null);
	}
}

/**************************************************************
 * InvalidRangeIndex
 **************************************************************/
class InvalidRangeIndex(TypeOfSize Row, TypeOfSize Column): Exception{
	import std.conv: toChars;

	/********************************************
	 * Constructors
	 ********************************************/
	@safe pure nothrow{
		shared static this(){
			foreach(elm; toChars(Row)) this.rowInStr ~= elm;
			this.strPostfix= "] is not allowed to access (" ~rowInStr;
			this.strPostfix ~= "x";
			foreach(elm; toChars(Column)) this.strPostfix ~= elm;
			this.strPostfix ~= ") matrix.";
		}

		this(in uint j){
			string buf= "Usage of index[1≤i≤" ~rowInStr ~", ";
			foreach(elm; toChars(j)) buf ~= elm;
			buf ~= strPostfix;
			super(buf);
		}

		this(in uint i, in uint j){
			string buf= "Usage of index[";
			foreach(elm; toChars(i)) buf ~= elm;
			buf ~= ", ";
			foreach(elm; toChars(j)) buf ~= elm;
			buf ~= strPostfix;
			super(buf, __FILE__, __LINE__, null);
		}
	}

private:
  static immutable string strPostfix;
	static immutable string rowInStr;
}

/**************************************************************
 * singular matrix
 **************************************************************/
class SingularMatrix: Exception{
	this() @safe pure nothrow @nogc{
		super("Singular matrix", __FILE__, __LINE__, null);
	}
}
