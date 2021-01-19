module sekitk.qvm.attribute;

import sekitk.qvm.common: TypeOfSize, MatrixType, MajorOrder;

struct MatrixAttribute(TypeOfSize Row, TypeOfSize Column, MatrixType Shape){
	import std.typecons: Ternary;
	import sekitk.qvm.common: ReturnTypeOfTranspose;

	@safe pure nothrow @nogc{
		/// copy constructor
		this(ref return scope inout typeof(this) other){}

		///
		this(inout scope bool invertible, inout scope bool orthogonal){
			this._invertible= invertible;
			this._orthogonal= orthogonal;
		}
	}

	@safe pure nothrow @nogc{
		MatrixAttribute!(Column, Row,
										 ReturnTypeOfTranspose!Shape) transpose() const{
		  auto result= typeof(return)();
			result._invertible= this._invertible;
			result._orthogonal= this._orthogonal;
			return result;
		}

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
		case MatrixType.permutation:
			result= Ternary.yes;
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
