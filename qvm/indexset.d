/******************************************************************************
 * Internal index of elements
 *
 * "ForwardRange"
 ******************************************************************************/
module sekitk.qvm.indexset;

import sekitk.qvm.common: TypeOfSize, TypeOfIndex, MatrixType, MajorOrder, matrixConstraint;

/**************************************************************
 * Mixin templates
 **************************************************************/
mixin template IndexSetBase(MatrixType Shape, IndexSetType IdxTyp)
if(Shape is MatrixType.zero ||
	 Shape is MatrixType.band3 ||
	 Shape is MatrixType.upperTri ||
	 Shape is MatrixType.lowerTri){
	import std.algorithm;
	import std.array;
	import std.range;
	import std.typecons: Tuple, tuple;
	import std.traits: Unqual, ParameterTypeTuple;
	import core.exception: onRangeError;

	/********************************************
	 * Constructors
	 ********************************************/
	this(in typeof(this) other) @safe pure nothrow @nogc{}

	/********************************************
	 * Requirement for a ForwardRange
	 ********************************************/
	@safe pure nothrow {
		TypeOfIndex front() @nogc const @property{
			onRangeError;
			return 0;
		}

		TypeOfIndex moveFront() @nogc{
			onRangeError;
			return 0;
		}

		void popFront() @nogc{onRangeError;}
		bool empty() @nogc const @property{return true;}
	}

	/********************************************
	 * Other methods
	 ********************************************/
	@safe pure nothrow @nogc const{
		TypeOfIndex back() const @property{
			onRangeError;
			return 0;
		}

		TypeOfIndex opIndex(in size_t index){
			onRangeError;
			return 0;
		}
		size_t length() @property{return 0;}
		alias opDollar= length;
	}

	@safe pure const{
		int opApply(scope int delegate(in TypeOfIndex) @safe pure nothrow dg){return 0;}
		int opApply(scope int delegate(in size_t, in TypeOfIndex) @safe pure nothrow dg){return 0;}
	}
}

/// ditto
mixin template IndexSetBase(TypeOfSize Row, TypeOfSize Column,
														MatrixType Shape,
														MajorOrder MatOdr,
														IndexSetType IdxTyp)
if(Shape is MatrixType.dense ||
	 Shape is MatrixType.band1 ||
	 Shape is MatrixType.band3 ||
	 Shape is MatrixType.upperTri ||
	 Shape is MatrixType.lowerTri){
	import core.exception: onRangeError;
	import std.algorithm: min;
	import std.traits: Unqual, ParameterTypeTuple;

	/********************************************
	 * Constructor
	 ********************************************/
	this(in typeof(this) other) @safe pure nothrow @nogc{
		this._idx= other._idx;
		this._itrStep= other._itrStep;
	}

	/********************************************
	 * Requirement for a InputRange
	 ********************************************/
	@safe pure nothrow{
		TypeOfIndex front() @nogc const @property{
			if(this.empty) onRangeError;
			return _idx;
		}

	  TypeOfIndex moveFront() @nogc{
		  mixin(INITIALIZE);
			return _idx;
		}

		void popFront() @nogc{
			if(empty) onRangeError;
			_idx= recurrRel(_idx, _itrStep++);
		}

		bool empty() @nogc const @property{
			return (_itrStep.totalStepRemain > 0)? false: true;
		}
	}

	/********************************************
	 * Other methods
	 ********************************************/
	@safe pure nothrow @nogc const{
		size_t length(){
			return _itrStep.totalStepRemain;
		}

		alias opDollar= length;
	}

	@safe{
		enum string IMPL= q{
			const size_t offset= _itrStep.totalStep;
			typeof(return) result= 0;
			TypeOfIndex idx= _idx;
			auto step= Unqual!(typeof(_itrStep))(_itrStep);
			foreach(_; offset..offset+_itrStep.totalStepRemain){
				result= dg(idx);
				idx= recurrRel(idx, step++);
			}
			return result;
		};

		int opApply(Dg)(scope Dg dg) if(ParameterTypeTuple!Dg.length == 1){
			mixin(IMPL);
		}
/+
		int opApply(DG: int delegate(in TypeOfIndex __applyArg0) @safe pure nothrow)(scope DG dg) pure nothrow const{
			mixin(IMPL);
		}

		int opApply(DG: int delegate(in T __applyArg0) @safe pure nothrow @nogc,
								T: TypeOfIndex)(scope DG dg) pure nothrow @nogc const{
			mixin(IMPL);
		}

		int opApply(DG: int delegate(ref T __applyArg0) @safe pure nothrow,
								T: TypeOfIndex)(scope DG dg) pure nothrow const{
			mixin(IMPL);
		}
+/
		int opApply(DG: int delegate(in size_t, in TypeOfIndex) @safe pure nothrow @nogc
								)(scope DG dg) pure nothrow @nogc const{
			const size_t offset= _itrStep.totalStep;
			TypeOfIndex idx= _idx;
			typeof(return) result= 0;
			auto step= Unqual!(typeof(_itrStep))(_itrStep);
			foreach(_; offset..offset+_itrStep.totalStepRemain){
				result= dg(step.totalStepCurr, idx);
				idx= recurrRel(idx, step++);
			}
			return result;
		}
	}

private:
	TypeOfIndex _idx= INDEX_INIT[0];
	IterationStep!(min(Row, Column), Shape, MatOdr, IdxTyp) _itrStep;
}

/**************************************************************
 * Index set for diagonal elements
 **************************************************************/
struct IndexSetDiag(TypeOfSize Row, TypeOfSize Column, MatrixType Shape, MajorOrder MatOdr)
if(Shape is MatrixType.zero){
	mixin IndexSetBase!(Shape, IndexSetType.diag);
}

/// ditto
struct IndexSetDiag(TypeOfSize Row, TypeOfSize Column, MatrixType Shape, MajorOrder MatOdr)
if(Shape !is MatrixType.zero &&
	 Shape !is MatrixType.permutation &&
	 matrixConstraint!(Row, Column, Shape)){
	import std.traits: Unqual;
	import std.algorithm: min;
	import std.range: iota, repeat;

	enum TypeOfIndex[1] INDEX_INIT= [0];
	enum string INITIALIZE= "_idx= INDEX_INIT[0];_itrStep.initialize;";

	mixin IndexSetBase!(Row, Column, Shape, MatOdr, IndexSetType.diag);

	/****************************
	 * Operators
	 ****************************/
	TypeOfIndex opIndex(in size_t itrStep) @safe pure nothrow @nogc const{
		import core.exception: onRangeError;
		import sekitk.integers.progression: sumFromZero;

		const i= itrStep+_itrStep.stepLocal;
	  size_t result;

		if(itrStep < length){
			final switch(Shape){
			case MatrixType.zero, MatrixType.permutation:
				assert(false);	// unreachable
			case MatrixType.dense:
				result= (MatOdr is MajorOrder.row)? i*(Column+1u) : i*(Row+1u);
				break;
			case MatrixType.band1:
				result= i;
				break;
			case MatrixType.band3:
				result= 3u*i;
				break;
			case MatrixType.upperTri:
				result= (MatOdr is MajorOrder.row)? i*Column-sumFromZero(i)+i : i+sumFromZero(i);
				break;
			case MatrixType.lowerTri:
				result= (MatOdr is MajorOrder.row)? i+sumFromZero(i) : i*Row-sumFromZero(i)+i;
			}
		}
		else{
			onRangeError();
		}
		return cast(typeof(return))(result);
	}

private:
	/********************************************
	 * $(I a)_(i+1)= f(a_i, i)
	 ********************************************/
	TypeOfIndex recurrRel(in TypeOfIndex idxCurr, in Unqual!(typeof(_itrStep)) itr) @safe pure nothrow @nogc const{
		typeof(return) result= idxCurr;

		final switch(Shape){
		case MatrixType.zero, MatrixType.permutation:
			assert(false);
			break;
		case MatrixType.dense:
			result += (MatOdr is MajorOrder.row)? Column+1u: Row+1u;
			break;
		case MatrixType.band1:
			result += 1u;
			break;
		case MatrixType.band3:
			result += 3u;
			break;
		case MatrixType.upperTri:
			result += (MatOdr is MajorOrder.row)? Row-itr.stepLocal: itr.stepLocal+2u;
/+
Matrix!(6, 6, Dense, Row)
 0  1  2  3  4  5
 *  6  7  8  9 10
 *  * 11 12 13 14
 *  *  * 15 16 17
 *  *  *  * 18 19
 *  *  *  *  * 20
idx= [0, 6, 11, 15, 18, 20]; d= 6, 5, 4, 3, 2

Matrix!(6, 6, Dense, Column)
 0  1  3  6 10 15
 *  2  4  7 11 16
 *  *  5  8 12 17
 *  *  *  9 13 18
 *  *  *  * 14 19
 *  *  *  *  * 20
idx= [0, 2, 5, 9, 14, 20]; d= 2, 3, 4, 5, 6
+/
			break;
		case MatrixType.lowerTri:
			result += (MatOdr is MajorOrder.row)? itr.stepLocal+2u: Column-itr.stepLocal;
		}
		return result;
	}
}

/**************************************************************
 * Index set for tridiagonal elements excluding diagonal ones
 **************************************************************/
struct IndexSetSubDiag(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr)
if(Shape is MatrixType.zero ||
	 Shape is MatrixType.band1 ||
	 Size < 3){
	mixin IndexSetBase!(Shape, IndexSetType.subDiag);
}

/// ditto
struct IndexSetSubDiag(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr)
if(Shape !is MatrixType.zero &&
	 Shape !is MatrixType.band1 &&
	 Size > 2u && matrixConstraint!(Size, Size, Shape)){
	import std.traits: Unqual;

	enum TypeOfIndex[_itrStep.MAX_ITR_SLICE] INDEX_INIT= (){
	  TypeOfIndex[_itrStep.MAX_ITR_SLICE] result;
		final switch(MatOdr){
		case MajorOrder.row:
			static if(Shape is MatrixType.dense) result= [1, Size];
			else static if(Shape is MatrixType.band3) result= [1, 2];
			else result= [1];
			break;
		case MajorOrder.column:
			static if(Shape is MatrixType.dense) result= [Size, 1];
			else static if(Shape is MatrixType.band3) result= [2, 1];
		  else result= [1];
		}
		return result;
	}();

	enum string INITIALIZE= "_idx= INDEX_INIT[0];_itrStep.initialize;";

	mixin IndexSetBase!(Size, Size, Shape, MatOdr, IndexSetType.subDiag);

	/****************************
	 * Operators
	 ****************************/
	TypeOfIndex opIndex(in size_t idx) @safe pure nothrow @nogc const{
		import sekitk.intop: partialSum;
		typeof(return) result;

		if(idx < length){
			auto pos= _itrStep.itrPosition(idx);
			result= INDEX_INIT[pos.stepSlice];
			if(pos.stepLocal > 0){
				final switch(Shape){
				case MatrixType.dense:
					result += (Size+1u)*pos.stepLocal;
					break;
				case MatrixType.band3:
					result += 3*pos.stepLocal;
					break;
				case MatrixType.upperTri, MatrixType.lowerTri:
					static if(Shape is MatrixType.upperTri && MatOdr is MajorOrder.row ||
										Shape is MatrixType.lowerTri && MatOdr is MajorOrder.column){
						result += pos.stepLocal*Size
							-sumFromZero(cast(TypeOfIndex)(pos.stepLocal-1));
					}
					else{
						result += summation(3u, pos.stepLocal+2u);
					}
					break;
				case MatrixType.zero, MatrixType.permutation, MatrixType.band1:
					assert(false);
				}
			}
		}
		else{
			onRangeError;
		}
		return cast(typeof(return))(result);
	}

private:
	TypeOfIndex recurrRel(in TypeOfIndex idxCurr, in Unqual!(typeof(_itrStep)) itr) @safe pure nothrow @nogc const{
		TypeOfIndex result;

		if(itr.isTurningPoint){
			result= INDEX_INIT[itr.stepSlice+1];
		}
		else{
			result= idxCurr;

			final switch(Shape){
			case MatrixType.dense:
				result += Size+1;
				break;
			case MatrixType.band3:
				result += 3u;
				break;
			case MatrixType.upperTri:
				result += (MatOdr is MajorOrder.row)? Size-itr.stepLocal: 3u+itr.stepLocal;
				break;
			case MatrixType.lowerTri:
				result += (MatOdr is MajorOrder.row)? 3u+itr.stepLocal: Size-itr.stepLocal;
				break;
			case MatrixType.zero, MatrixType.permutation, MatrixType.band1:
				assert(false);
			}
		}
		return result;
	}
}

/**************************************************************
 * Index set for strict upper triangular elements
 **************************************************************/
struct IndexSetStrictTriR(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr)
if(Shape is MatrixType.zero ||
	 Shape is MatrixType.band1 ||
	 Shape is MatrixType.lowerTri ||
	 Size < 2 || !matrixConstraint!(Size, Size, Shape)){
	mixin IndexSetBase!(Shape, IndexSetType.strictUpperTri);
}

/// ditto
struct IndexSetStrictTriR(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr)
if(Shape !is MatrixType.zero &&
	 Shape !is MatrixType.band1 &&
	 Shape !is MatrixType.lowerTri &&
	 Size > 1 && matrixConstraint!(Size, Size, Shape)){
	import std.traits: Unqual;
	import std.range: iota, recurrence, repeat;
	import sekitk.integers.progression: sumFromZero;

	/********************************************
	 * Manifest constants
	 ********************************************/
	static if(Shape is MatrixType.dense || Shape is MatrixType.upperTri){
		enum TypeOfIndex[_itrStep.MAX_ITR_SLICE] INDEX_INIT= (){
			import std.array: staticArray;
			enum TypeOfSize LEN= _itrStep.MAX_ITR_SLICE;
			TypeOfIndex[LEN] num;
			static if(MatOdr is MajorOrder.row){
				num= recurrence!((a, n) => cast(TypeOfIndex)(a[n-1]+Size+1u)
												 )(TypeOfIndex(1u)).staticArray!LEN;
				static if(Shape is MatrixType.upperTri){
					foreach(i; 1u..LEN) num[i] -= sumFromZero(i);
				}
			}
			else static if(MatOdr is MajorOrder.column){
				static if(Shape is MatrixType.dense){
					num= recurrence!((a, n) => a[n-1]+Size
													 )(TypeOfIndex(Size)).staticArray!LEN;
				}
				else{
					num= recurrence!((a, n) => cast(TypeOfIndex)(sumFromZero(n+1))
													 )(TypeOfIndex(1u)).staticArray!LEN;
				}
			}
			else{
				num= recurrence!((a, n) => a[i-1]+2*(Size-n)
												 )(TypeOfIndex(Size)).staticArray!LEN;
				/+num[0]= Size;
				num[i]= num[i-1]+2*(Size-i);+/
			}
			return num;
		}();
	}
	else{	// Band3
		static if(MatOdr is MajorOrder.row){
			enum TypeOfIndex[1] INDEX_INIT= 1;
		}
		else static if(MatOdr is MajorOrder.column){
			enum TypeOfIndex[1] INDEX_INIT= 2;
		}
		else{
			enum TypeOfIndex[1] INDEX_INIT= Size;
		}
	}

	enum string INITIALIZE= "_idx= INDEX_INIT[0];_itrStep.initialize;";

	mixin IndexSetBase!(Size, Size, Shape, MatOdr, IndexSetType.strictUpperTri);

	/****************************
	 * Operators
	 ****************************/
	TypeOfIndex opIndex(in size_t idx) @safe pure nothrow @nogc const{
		typeof(return) result;
		if(idx < length){
			auto pos= _itrStep.itrPosition(idx);
			result= INDEX_INIT[pos.stepSlice];
			if(pos.stepLocal > 0){
				final switch(Shape){
				case MatrixType.band3:
					break;
				case MatrixType.dense, MatrixType.upperTri:
					result += pos.stepLocal;
					break;
				case MatrixType.zero, MatrixType.permutation,
					MatrixType.band1, MatrixType.lowerTri:
					assert(false);
				}
			}
		}
		else{
			onRangeError;
		}
		return result;
	}

private:
	TypeOfIndex recurrRel(in TypeOfIndex idxCurr, in Unqual!(typeof(_itrStep)) itr) @safe pure nothrow @nogc const{
		typeof(return) result;

		if(itr.isTurningPoint){
			result= INDEX_INIT[itr.stepSlice+1];
		}
		else{
			result= idxCurr;

			final switch(Shape){
			case MatrixType.dense, MatrixType.upperTri, MatrixType.band3:
				result += 1u;
				break;
			case MatrixType.zero, MatrixType.permutation, MatrixType.band1, MatrixType.lowerTri:
				assert(false);
			}
		}

		return result;
	}
}

/**************************************************************
 * strict lower triangular elements
 **************************************************************/
struct IndexSetStrictTriL(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr)
if(Shape is MatrixType.zero ||
	 Shape is MatrixType.band1 ||
	 Shape is MatrixType.upperTri ||
	 Size < 2){
	mixin IndexSetBase!(Shape, IndexSetType.strictLowerTri);
}

/// ditto
struct IndexSetStrictTriL(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr)
if(Shape !is MatrixType.zero &&
	 Shape !is MatrixType.band1 &&
	 Shape !is MatrixType.upperTri &&
	 Size > 1 && matrixConstraint!(Size, Size, Shape)){
	import std.traits: Unqual;
	import std.range: iota, repeat;
	import sekitk.integers.progression: sumFromZero;

	/********************************************
	 * Manifest constants
	 ********************************************/
	static if(Shape is MatrixType.dense || Shape is MatrixType.lowerTri){
		enum TypeOfIndex[_itrStep.MAX_ITR_SLICE] INDEX_INIT= (){
			TypeOfIndex[_itrStep.MAX_ITR_SLICE] num;
			static if(MatOdr is MajorOrder.column){
				num[0]= 1;
				static foreach(i; 1u..num.length){
					num[i]= cast(TypeOfIndex)(num[i-1]+(Size+1u));
					static if(Shape is MatrixType.lowerTri) num[i] -= sumFromZero(i);
				}
			}
			else static if(MatOdr is MajorOrder.row){
				static if(Shape is MatrixType.dense){
					num[0]= Size;
					static foreach(TypeOfIndex i; 1..num.length) num[i]= cast(TypeOfIndex)(num[i-1]+Size);
				}
				else{
					static foreach(TypeOfIndex i; 0..num.length) num[i]= sumFromZero(cast(TypeOfIndex)(i+1u));
				}
			}
			else{
				num[0]= Size;
				static foreach(TypeOfIndex i; 1..num.length) num[i]= num[i-1]+2*(Size-i);
			}
			return num;
		}();
	}
	else{	// Band3
		static if(MatOdr is MajorOrder.column){
			enum TypeOfIndex[1] INDEX_INIT= 1;
		}
		else static if(MatOdr is MajorOrder.row){
			enum TypeOfIndex[1] INDEX_INIT= 2;
		}
		else{
			enum TypeOfIndex[1] INDEX_INIT= Size;
		}
	}

	enum string INITIALIZE= "_idx= INDEX_INIT[0];_itrStep.initialize;";

	mixin IndexSetBase!(Size, Size, Shape, MatOdr, IndexSetType.strictLowerTri);

	/********************************************
	 * Other methods
	 ********************************************/
	TypeOfIndex opIndex(in size_t idx) @safe pure nothrow @nogc const{
		typeof(return) result;
		if(idx < length){
			auto pos= _itrStep.itrPosition(idx);
			result= INDEX_INIT[pos.stepSlice];
			if(pos.stepLocal > 0){
				final switch(Shape){
				case MatrixType.band3:
					break;
				case MatrixType.dense, MatrixType.upperTri:
					result += pos.stepLocal;
					break;
				case MatrixType.zero, MatrixType.permutation, MatrixType.band1, MatrixType.lowerTri:
					assert(false);
				}
			}
		}
		else{
			onRangeError;
		}
		return result;
	}

private:
	TypeOfIndex recurrRel(in TypeOfIndex idxCurr, in Unqual!(typeof(_itrStep)) itr) @safe pure nothrow @nogc const{
		typeof(return) result;

		if(itr.isTurningPoint){
			result= INDEX_INIT[itr.stepSlice+1];
		}
		else{
			result= idxCurr;

			final switch(Shape){
			case MatrixType.dense, MatrixType.lowerTri, MatrixType.band3:
				result += 1u;
				break;
			case MatrixType.zero, MatrixType.permutation, MatrixType.band1, MatrixType.upperTri:
				assert(false);
			}
		}

		return result;
	}
}

/**********************************************
 * transpose
 **********************************************/
struct IndexSetTranspose(TypeOfSize Row, TypeOfSize Column, MatrixType Shape, MajorOrder MatOdr)
if(Shape !is MatrixType.zero && matrixConstraint!(Row, Column, Shape)){
	import core.exception: onRangeError;
	import std.algorithm: sum;
	import std.array: staticArray;
	import std.range: iota, recurrence, repeat, retro;

	enum TypeOfSize MAX_ITR_SLICE= (MatOdr is MajorOrder.row)? Row: Column;
	enum TypeOfSize[MAX_ITR_SLICE] MAX_ITR_LOCAL= (){
		static if(Shape is MatrixType.dense){
			static if(MatOdr is MajorOrder.row){
				auto temp= Column.repeat(MAX_ITR_SLICE);
			}
			else static if(MatOdr is MajorOrder.column){
				auto temp= Row.repeat(MAX_ITR_SLICE);
			}
		}
		else static if(Shape is MatrixType.band1){
			auto temp= 1.repeat(MAX_ITR_SLICE);
		}
		else static if(Shape is MatrixType.band3){
		  TypeOfSize[] temp= new TypeOfSize[MAX_ITR_SLICE];
			temp[0]= 2;
			temp[1..MAX_ITR_SLICE-1]= 3;
			temp[MAX_ITR_SLICE-1]= 2;
		}
		else static if(Shape is MatrixType.upperTri ||
									 Shape is MatrixType.lowerTri){
			auto temp= iota!(TypeOfSize, TypeOfSize)(1u, Size+1u).retro;
		}
		else assert(false);
		return temp.staticArray!MAX_ITR_SLICE;
	}();
	enum TypeOfIndex MAX_ITR_TOTAL= MAX_ITR_LOCAL[].sum;

	/****************************
	 * InputRange
	 ****************************/
	@safe pure nothrow @nogc{
		TypeOfIndex front() const @property{
			return idx;
		}

		void popFront(){
			if(stepLocal == MAX_ITR_LOCAL[stepSlice]-1){
				++stepSlice;
				stepLocal= 0;
			}
			else{
				++stepLocal;
			  if(stepTotal >= MAX_ITR_TOTAL) onRangeError;
			}
			++stepTotal;
			idx= recurrRel(idx, stepSlice, stepLocal);
		}

		bool empty() const @property{
			if(stepTotal < MAX_ITR_TOTAL) return false;
			else return true;
		}
	}

	size_t length() @safe pure nothrow @nogc const @property{
		return MAX_ITR_TOTAL-stepTotal;
	}

	alias opDollar= length;

	auto seed(){return idxBeg;}
private:
	TypeOfIndex recurrRel(in TypeOfIndex _idx, in TypeOfSize _stepSlice, in TypeOfSize _stepLocal) @safe pure nothrow @nogc const{
		typeof(return) result;

		if(_stepLocal == 0 && _stepSlice < MAX_ITR_SLICE) result= idxBeg[_stepSlice];
		else{
			result= _idx;
			final switch(Shape){
			case MatrixType.dense:
				final switch(MatOdr){
				case MajorOrder.row:
					result += Row;
					break;
				case MajorOrder.column:
					result += Column;
				}
				break;
			case MatrixType.band3:
				result += 2;
				break;
			case MatrixType.upperTri, MatrixType.lowerTri:
				result += _stepSlice+_stepLocal;
				break;
			case MatrixType.zero, MatrixType.permutation, MatrixType.band1:
				assert(false);
			}
		}
		return result;
	}

	TypeOfSize stepSlice= 0u;
	TypeOfSize stepLocal= 0u;
	TypeOfIndex stepTotal= 0u;
	TypeOfIndex idx= 0u;
	static if(Shape is MatrixType.dense){
		static if(MatOdr is MajorOrder.row){
			TypeOfIndex[Row] idxBeg= iota(Row).staticArray!Row;
		}
		else static if(MatOdr is MajorOrder.column){
			TypeOfIndex[Column] idxBeg= iota(Column).staticArray!Column;
		}
	}
	else{
		enum TypeOfSize Size= Row;
		TypeOfIndex[Size] idxBeg= (){
			import std.array: array, staticArray;
			import std.range: take;

			static if(Shape is MatrixType.band1){
				auto idxSet= iota(0u, Size);
			}
			else static if(Shape is MatrixType.band3){
			  TypeOfIndex[] idxSet= new TypeOfIndex[Size];
				idxSet[0]= 0u;
				idxSet[1..Size]= iota!(TypeOfIndex, TypeOfIndex, TypeOfIndex)(1, 3*Size-4, 3).array;
			}
			else{
				auto idxSet= recurrence!"a[n-1]+n+1"(0, 2).take(Size);
			}
			return idxSet.staticArray!Size;
		}();
	}
}

/**************************************************************
 * Private objects
 **************************************************************/
private enum IndexSetType{diag, subDiag, strictUpperTri, strictLowerTri}

private struct IterationStep(TypeOfSize Size, MatrixType Shape, MajorOrder MatOdr, IndexSetType Typ)
if((Size > 0 && Typ is IndexSetType.diag) ||
	 (Size > 2 && Typ is IndexSetType.subDiag) ||
	 (Size > 1 && Typ is IndexSetType.strictUpperTri) ||
	 (Size > 1 && Typ is IndexSetType.strictLowerTri)){
	import core.exception: onRangeError;
	import std.algorithm: sum;
	import std.range: iota, retro;

	/********************************************
	 * Manifest constants
	 ********************************************/
	static if(Typ is IndexSetType.diag){
		enum TypeOfSize MAX_ITR_SLICE= 1;
		enum TypeOfSize[MAX_ITR_SLICE] MAX_ITR_LOCAL= [Size];
	}
	else static if(Typ is IndexSetType.subDiag){
		enum TypeOfSize MAX_ITR_SLICE= (){
			TypeOfSize result;
			final switch(Shape){
			case MatrixType.dense, MatrixType.band3:
				result= 2u;
				break;
			case MatrixType.upperTri, MatrixType.lowerTri:
				result= 1u;
				break;
			case MatrixType.zero, MatrixType.permutation, MatrixType.band1:
				assert(false);
			}
			return result;
		}();
		enum TypeOfSize[MAX_ITR_SLICE] MAX_ITR_LOCAL= Size-1;
	}
	else static if(Typ is IndexSetType.strictUpperTri ||
								 Typ is IndexSetType.strictLowerTri){
		enum TypeOfSize MAX_ITR_SLICE= Size-1;
		enum TypeOfSize[MAX_ITR_SLICE] MAX_ITR_LOCAL= (){
			import std.array: staticArray;
			static if((Typ is IndexSetType.strictUpperTri && MatOdr is MajorOrder.row) ||
								(Typ is IndexSetType.strictLowerTri && MatOdr is MajorOrder.column)){
				auto temp= iota(1u, Size).retro;
			}
			else{
				auto temp= iota(1u, Size);
			}
			return temp.staticArray!MAX_ITR_SLICE;
		}();
	}

	enum TypeOfSize MAX_TOTAL_STEP= MAX_ITR_LOCAL[].sum;

	/********************************************
	 * Constructor
	 ********************************************/
	this(in typeof(this) other) @safe pure nothrow @nogc{
		this.stepSlice= other.stepSlice;
		this.stepLocal= other.stepLocal;
	}

	@safe pure nothrow @nogc{
		/// Operator
		typeof(this) opUnary(string Inc)() if(Inc == "++"){
			import core.exception: onRangeError;
			if(this.isTurningPoint){
				if(stepSlice == MAX_ITR_SLICE) onRangeError;
				else{
					++ stepSlice;
					stepLocal= 0;
				}
			}
			else ++ stepLocal;

			return this;
		}

		void initialize(){
			this.stepSlice= 0u;
			this.stepLocal= 0u;
		}
	}

	@safe pure nothrow @nogc const{
		immutable(TypeOfSize) totalSlice() @property{
			return MAX_ITR_SLICE;
		}

		immutable(TypeOfIndex) totalStep() @property{
			TypeOfIndex num= 0u;
			foreach(slc; 0..MAX_ITR_SLICE) num += MAX_ITR_LOCAL[slc];
			return num;
		}

		auto itrPosition(in size_t totalStep){
			import std.algorithm: maxElement;
			import std.typecons: tuple;

			if(totalStep >= this.totalStep) onRangeError;

			TypeOfSize posSlice= MAX_ITR_SLICE;
			TypeOfSize posLocal= cast(TypeOfSize)totalStep;	// temp
			TypeOfIndex numBef;
			TypeOfIndex numCurr= 0u;
			foreach(TypeOfSize slc; 0..MAX_ITR_SLICE){
				numBef= numCurr;
				numCurr += MAX_ITR_LOCAL[slc];
				if(totalStep < numCurr){
					posSlice= slc;
					posLocal -= numBef;
					break;
				}
				else continue;
			}
			assert(posSlice < MAX_ITR_SLICE);
			return tuple!("stepSlice", "stepLocal")(posSlice, posLocal);
		}

		bool isTurningPoint() @property{
			bool result;
			static if(MAX_ITR_SLICE == 1) result= false;
			else result= (stepSlice < MAX_ITR_SLICE-1 &&
							stepLocal == MAX_ITR_LOCAL[stepSlice]-1)? true: false;
			return result;
		}

		TypeOfSize totalStepCurr() @property{
			typeof(return) result= 0u;
			foreach(idx; 0u..stepSlice) result += MAX_ITR_LOCAL[idx];
		  result += stepLocal;

			return result;
		}

		TypeOfSize totalStepRemain() @property{
			return cast(TypeOfSize)(MAX_TOTAL_STEP-this.totalStepCurr());
		}
	}

	TypeOfSize stepSlice;
	TypeOfSize stepLocal;
}
