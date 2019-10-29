module sekitk.attribute;

import sekitk.base: TypeOfSize, MatrixType, MajorOrder;

struct MatrixAttribute(TypeOfSize Row, TypeOfSize Column, MatrixType Shape){
	import std.typecons: Ternary;

	this(inout scope bool invertible, inout scope bool orthogonal) @safe pure nothrow @nogc{
		this._invertible= invertible;
		this._orthogonal= _orthogonal;
	}

	@safe pure nothrow @nogc{
		bool isInvertible() const{
			return (_invertible is Ternary.yes);
		}

		bool isSingular() const{
			return (_invertible is Ternary.no);
		}

		void setInvertibility(inout scope bool state){
			_invertible= state? Ternary.yes: Ternary.no;
		}
	}

private:
	Ternary _invertible= (){
		Ternary result;
		final switch(Shape){
		case MatrixType.zero:
			result= Ternary.no;
			break;
		case MatrixType.dense,
			MatrixType.band1,
			MatrixType.band3,
			MatrixType.upperTri,
			MatrixType.lowerTri:
			result= (Row == Column)? Ternary.unknown: Ternary.no;
		}
		return result;
	}();

	Ternary _orthogonal;//= _invertible;
}
