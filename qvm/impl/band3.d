/*****************************************************************************
 * Tridiagonal matrix
 *****************************************************************************/
module sekitk.qvm.impl.band3;

import sekitk.qvm.base: TypeOfSize, MajorOrder;

/*************************************************************
 * TriDiagonalMatrix
 *************************************************************/
mixin template Band3(T, real Threshold,
		     TypeOfSize Size, MajorOrder MatOdr)
if(Size > 2u){
  import std.array: array, staticArray;
  import std.range: iota, take, retro, zip;
  import mathematic.basic: Identity;
  import sekitk.qvm.base;
  import sekitk.qvm.indexset;

public:
  /*******************************************
   * Operators
   *******************************************/
  @safe pure nothrow const{
    // add & sub (different shape)
    // band3 pm band1
    Matrix!(Row, Column,
	    MatrixType.band3) opBinary(string Op,
				       MatrixType ShapeR: MatrixType.band1
				       )(in Matrix!(Row, Column, ShapeR) rhs)
    if(isPlusOrMinusSign!Op){
      enum Size_t LEN= arrayLength!(Row, Column, MatrixType.band3);
      T[LEN] num= this._values[];
      auto idxSet= IndexSetDiag!TmplArgsOfThis();

      foreach(scope idxSrc, scope idxDest;
	      zip(iota(LEN), idxSet)) mixin("num[idxDest] " ~Op ~"= rhs.v[idxSrc];");
      return new typeof(return)(num);
    }

    // band3 pm dense, band3 pm upper triangular, band3 pm lower triangular
    Matrix!(Row, Column,
	    MatrixType.dense) opBinary(string Op,
				       MatrixType ShapeR)(in Matrix!(Row, Column, ShapeR) rhs)
    if(isPlusOrMinusSign!PM && (ShapeR is MatrixType.dense ||
				ShapeR is MatrixType.upperTri ||
				ShapeR is MatrixType.lowerTri)){
      enum MatrixType ShapeDest= MatrixType.dense;
      enum TypeOfIndex LEN= arrayLength!(Row, Column, MatrixType.dense);
      T[LEN] num= VALUE_ZERO;
      auto idxSetDest0= IndexSetDiag!(Row, Column, ShapeDest, MatOdr)();
      auto idxSetDest1= IndexSetSubDiag!(Size, ShapeDest, MatOdr)();
      auto idxSetSrc= arrayLength!(Row, Column, Shape).iota();
      foreach(scope idxSrc, scope idxDest;
	      zip(idxSetSrc, merge(idxSet0, idxSet1))) num[idxDest]= this._values[idxSrc];

      static if(ShapeR is MatrixType.dense){
	foreach(scope idx; 0u..LEN) mixin("num[idx]" ~Op ~"rhs.v[idx];");
      }
      else{
	auto idxSetDest0= IndexSetDiag!(Row, Column, ShapeDest, MatOdr)();
	static if(ShapeR is MatrixType.upperTri){
	  auto idxSetDest1= IndexSetStrictTriR!(Size, ShapeDest, MatOdr)();
	}
	else{
	  auto idxSetDest1= IndexSetStrictTriL!(Size, ShapeDest, MatOdr)();
	}
	auto idxSetSrc= arrayLength!(Row, Column, ShapeR).iota();
	foreach(scope idxSrc, scope idxDest;
		zip(idxSetSrc, merge(idxSetDest0, idxSetDest1)) )
	  mixin("num[idxDest]" ~Op ~"rhs.v[idxSrc];");
      }

      return new typeof(return)(temp);
    }

		
    /*****************************************
     * Matrix multiplication
     * band3.opBinary!"*"(band1)
     *****************************************/
    Matrix!(Row, ColumnR,
	    MatrixType.band3) opBinary(string Op: "*", TypeOfSize ColumnR,
				       MatrixType ShapeR: MatrixType.band1
				       )(in Matrix!(Column, ColumnR, ShapeR) rhs){
      enum TypeOfIndex LEN= arrayLength!(Row, ColumnR, MatrixType.band3);
      T[LEN] num;
      size_t st= 0, en;
      TypeOfIndex[] idxSetDest= iota!TypeOfIndex(LEN).array;
      TypeOfIndex[] idxSetRhs= (){
	auto seed= iota!TypeOfIndex(Size);
	TypeOfIndex[] result;

	result= seed.take(2).array;
	foreach(scope i; 1u..Size-1u){
	  en= st+3;
	  result ~= seed[st..en].array;
	  ++st;
	}
	result ~= seed[$-2..$].array;
	auto temp= result;

	return result;
      }();

      foreach(scope idxRhs, scope idxDest;
	      zip(idxSetRhs, idxSetDest)) num[idxDest]=
					    this._values[idxDest]*rhs._values[idxRhs];

      return new typeof(return)(num);
    }

    /*****************************************
     * Matrix multiplication
     *
     *****************************************/
    Matrix!(Row, ColumnR,
	    MatrixType.dense) opBinary(string Op: "*", TypeOfSize ColumnR,
				       MatrixType ShapeR
				       )(in Matrix!(Column, ColumnR, ShapeR) rhs)
    if(ShapeR is MatrixType.dense ||
       ShapeR is MatrixType.band3 ||
       ShapeR is MatrixType.upperTri ||
       ShapeR is MatrixType.lowerTri){
      import mathematic.progression: sumFromZero;
      alias idxMapDense= indexMap!(Row, ColumnR, MatrixType.dense, MatOdr);
      alias idxMapBand3= indexMap!TmplArgsOfThis;
      enum TypeOfIndex LEN= arrayLength!(Row, ColumnR, MatrixType.dense);
      T[LEN] num= void;

      final switch(ShapeR){
      case MatrixType.dense:	// tridiagonal * dense
	auto idxSetDest= iota!TypeOfIndex(Row*ColumnR);

	foreach(scope j; iota!TypeOfSize(ColumnR)){
	  num[idxSetDest.front]= this._values[0]*rhs._values[j]
	    +this._values[1]*rhs._values[j+ColumnR];
	  idxSetDest.popFront;
	}
	foreach(scope i; 1u..Row-1u){
	  foreach(scope j; 0u..ColumnR){
	    num[idxSetDest.front]= this._values[2+3*(i-1)]*rhs._values[j+(i-1)*ColumnR]
	      +this._values[3+3*(i-1)]*rhs._values[j+i*ColumnR]
	      +this._values[4+3*(i-1)]*rhs._values[j+(i+1)*ColumnR];
	    idxSetDest.popFront;
	  }
	}
	foreach(scope j; ColumnR.iota){
	  num[idxSetDest.front]= this._values[2+3*(Row-2)]*rhs._values[j+(Row-2)*ColumnR]
	    +this._values[3+3*(Row-1 -1)]*rhs._values[j+(Row-1)*ColumnR];
	  idxSetDest.popFront;
	}
	break;
      case MatrixType.band3:	// tridiagonal * tridiagonal
	{
	  TypeOfIndex[] idxSet;
	  TypeOfIndex st;
	  idxSet= iota!(TypeOfIndex, TypeOfIndex)(3u, Size).array;
	  static if(Size > 4){
	    st= 1u;
	    foreach(_; 0u..Size-3u){
	      st += Size+1u;
	      idxSet ~= iota!TypeOfIndex(st, st+Size-4u).array;//st= 3+Size-1, st..st+Size-4
	    }
	  }
	  st= Size*(Size-1u);
	  idxSet ~= iota!TypeOfIndex(st, cast(TypeOfIndex)(st+Size-3u)).array;
	  assert(idxSet.length == 2*sumFromZero(Size-3u));
	}
	// corner elements
	foreach(idx; idxSet) num[idx]= VALUE_ZERO;
	/+
size=4: total= 2*1
3,
12= Size*(Size-1)

size=5: total= 6= 2+2*2
3, 4,
9,
15,
20, 21

size=6: total= 12= 2+2*2+2*3
3, 4, 5,
10, 11,
17, 18,
24, 25,
30, 31, 32

size=7: total= 20= 2+4+6+8
3, 4, 5, 6,
11, 12, 13,
19, 20, 21,
27, 28, 29,
35, 36, 37,
42, 43, 44, 45
+/
	foreach(scope TypeOfSize j; 0u..Size-2u){
	  num[idxMapDense(j, j+2u)]=
	    this._values[idxMapBand3(j, j+1u)]*rhs._values[idxMapBand3(j+1u, j+2u)];
	  num[idxMapDense(j+2u, j)]=
	    this._values[idxMapBand3(j+2u, j+1u)]*rhs._values[idxMapBand3(j+1u, j)];
	}
	foreach(scope TypeOfSize j; 0u..Size-1u){
	  num[idxMapDense(j, j+1u)]=
	    this._values[idxMapBand3(j, j)]*rhs._values[idxMapBand3(j, j+1u)]
	    +this._values[idxMapBand3(j, j+1u)]*rhs._values[idxMapBand3(j+1u, j+1u)];
	  num[idxMapDense(j+1u, j)]=
	    this._values[idxMapBand3(j+1u, j)]*rhs._values[idxMapBand3(j, j)]
	    +this._values[idxMapBand3(j+1u, j+1u)]*rhs._values[idxMapBand3(j+1u, j)];
	}
	num[idxMapDense(0u, 0u)]= this._values[0]*rhs._values[0]
	  +this._values[idxMapBand3(0u, 1u)]*rhs._values[idxMapBand3(1u, 0u)];
	num[idxMapDense(Size-1u, Size-1u)]=
	  this._values[idxMapBand3(Size-1u, Size-2u)]*rhs._values[idxMapBand3(Size-2u, Size-1u)]
	  +this._values[idxMapBand3(Size-1u, Size-1u)]*rhs._values[idxMapBand3(Size-1u, Size-1u)];
	foreach(scope TypeOfSize j; 1u..Size-1u){
	  num[idxMapDense(j, j)]=
	    this._values[idxMapBand3(j, j-1u)]*rhs._values[idxMapBand3(j-1u, j)]
	    +this._values[idxMapBand3(j, j)]*rhs._values[idxMapBand3(j, j)]
	    +this._values[idxMapBand3(j, j+1u)]*rhs._values[idxMapBand3(j+1u, j)];
	}
	break;
      case MatrixType.upperTri:	// tridiagonal * upper triangular
	break;
      case MatrixType.lowerTri:	// tridiagonal * lower triangular
	break;
      case MatrixType.zero, MatrixType.permutation, MatrixType.band1:
	assert(false);
      }
      return new typeof(return)(num);
    }
  }

private:
  @safe pure nothrow @nogc const{
    /*****************************************
     * arrayDense
     *
     * This function generates the internal array of Matrix!(Size, Size, MatrixType.Dense, MatOdr)
     *****************************************/
    T[arrayLength!(Size, Size, MatrixType.dense)] arrayDense(){
      import std.array: array;
      import std.algorithm: merge;

      typeof(return) num= VALUE_ZERO;
      size_t j= 0;
      auto idxSet0= IndexSetDiag!(Size, Size, MatrixType.dense, MatOdr)();
      auto idxSet1= IndexSetSubDiag!(Size, MatrixType.dense, MatOdr)();
      foreach(scope TypeOfIndex i; merge(idxSet0, idxSet1)) num[i]= this._values[j++];

      return num;
    }
 
    /*****************************************
     * Determinant
     *
     * O(n)
     * f_n= a_n*f_{n-1}-c_{n-1}*b_{n-1}*f*_{n-2}
     *
     * See_Also:
     * https://en.wikipedia.org/wiki/Tridiagonal_matrix
     *****************************************/
    T detImpl(){
      typeof(return) result= void;
      T f_m1= Identity!(T, "*"), f_m2= Identity!(T, "+");
      auto idxSet= IndexSetDiag!(Row, Column, Shape, MatOdr)();

      foreach(scope TypeOfSize i; 1u..Size){
	result= idxSet[i+1u]*f_m1
	  -_values[idxLowerElm(i)]*_values[idxUpperElm(i)]*f_m2;
	f_m2= f_m1;
	f_m1= result;
      }
      return result;
    }

    /*****************************************
     * Inverse matrix
     *****************************************/
    T[arrayLength!(Size, Size, MatrixType.band3)] inverseImpl(){
      alias idxMap= indexMap!TemplateArgs;
      immutable TypeOfIndex[Size] idxSetA= (){
        auto temp= IndexSetDiag!TemplateArgs();
	return temp.staticArray!Size;
      }();
      immutable TypeOfIndex[Size-1] idxSetB= (){
        auto temp= IndexSetStrictTriR!(Size, Shape, MatOdr)();
	return temp.staticArray!(Size-1);
      }();
      immutable TypeOfIndex[Size-1] idxSetC= (){
        auto temp= IndexSetStrictTriL!(Size, Shape, MatOdr)();
	return temp.staticArray!(Size-1);
      }();

      const T[Size+2] phi= (){
        T[Size+2] temp;
	temp[Size]= Identity!(T, "*");
	temp[Size-1]= this._values[idxSetA[$-1u]];
	foreach_reverse(scope i; 1u..Size-1u) temp[i]= this._values[idxSetA[i]]*temp[i+1];
	return temp;
      }();

      const T[Size+1] theta= (){
        T[Size+1] temp= void;
	temp[0]= Identity!(T, "*");
	temp[1]= this._values[idxMap(0u, 0u)];
	foreach(scope i; 2u..Size+1u)
	  temp[i]= this._values[idxSetA[i]]*temp[i-1]
	    -this._values[idxSetB[i-1]]*this._values[idxSetC[i-1]]*temp[i-2];
	return temp;
      }();

      typeof(return) num= void;

      foreach(scope i; 0..Size){
	foreach(scope j; 0..Size){
	  if(i <= j){
	    num[idxMap(i, j)]= theta[i-1]*phi[j+1];
	    foreach(scope k; i..j-1u) num[idxMap(i, j)] *= this._values[idxSetB[k]];
	  } 
	  else{
	    num[idxMap(i, j)]= theta[j-1]*phi[i+1];
	    foreach(scope k; j..i-1u) num[idxMap(i, j)] *= this._values[idxSetC[k]];
	  }
	  num[idxMap(i, j)] /= theta[$-1];
	}
      }
      return num;
    }
  }

  /*******************************************
   * sub diagonal elements
   *******************************************/
  @safe pure nothrow @nogc const{
    TypeOfIndex idxUpperElm(in TypeOfSize index)
    in(index > 0 && index < Size){
      return cast(TypeOfIndex)(3u*index -(MatOdr is MajorOrder.row)? 2u: 1u);
    }

    TypeOfIndex idxLowerElm(in TypeOfSize index)
    in(index > 0 && index < Size){
      return cast(TypeOfIndex)(3u*index -(MatOdr is MajorOrder.row)? 1u: 2u);
    }
  }
}
