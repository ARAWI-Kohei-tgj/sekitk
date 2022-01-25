/******************************************************************************
 * dense matrix
 ******************************************************************************/
module sekitk.qvm.impl.dense;

import sekitk.qvm.common: TypeOfSize, MajorOrder;

/**************************************************************
 * common methods of dense matrix both rectangular and square
 **************************************************************/
/+
mixin template DenseCommon(T, TypeOfSize Row, TypeOfSize Column, MajorOrder MatOdr)
if(Row > 0 && Column > 0){
+/
enum string DENSE_COMMON= q{
  import std.traits: Unqual, ParameterTypeTuple;
  import std.array: array;
  import std.range: repeat, iota, zip;
  import sekitk.base: Identity;

  import sekitk.qvm.common;
  import sekitk.qvm.indexset;

private @safe:
  /********************************************
   * Operators
   ********************************************/
  pure nothrow const{
    // add & sub (different shape)
    TypeOfThis opBinary(string Op,
			MatrixType ShapeR: MatrixType.band1)(in Matrix!(Row, Column,
									ShapeR, MatOdr) rhs)
    if(isPlusOrMinusSign!Op){
      TypeOfInternalArray num= this._values[];
      TypeOfIndex j= 0u;
      auto idxSet= IndexSetDiag!TemplateArgs();
      foreach(scope i; idxSet) mixin("num[i] " ~PM ~"= rhs.v[j++];");

      return new typeof(return)(num);
    }

    // other shape
    TypeOfThis opBinary(string Op,											MatrixType ShapeR)(in Matrix!(Row, Column, ShapeR, MatOdr) rhs)
    if(isPlusOrMinusSign!Op && ShapeR != Shape){
      import std.algorithm: merge;
      TypeOfInternalArray num= this._values[];
      TypeOfIndex j= 0u;

      auto idxSet0= IndexSetDiag!TemplateArgs();

      static if(ShapeR is MatrixType.band3){
	auto idxSet1= IndexSetSubDiag!(Size, Shape, MatOdr)();
	foreach(idx; merge(idxSet0, idxSet1)) mixin("num[idx] " ~PM ~"= rhs._values[j++];");
      }
      else static if(ShapeR is MatrixType.upperTri){
	auto idxSet1= IndexSetStrictTriR!(Size, Shape, MatOdr)();
	foreach(idx; merge(idxSet0, idxSet1)) mixin("num[idx] " ~PM ~"= rhs._values[j++];");
      }
      else static if(ShapeR is MatrixType.lowerTri){
	auto idxSet1= IndexSetStrictTriL!(Size, Shape, MatOdr)();
	foreach(idx; merge(idxSet0, idxSet1)) mixin("num[idx] " ~PM ~"= rhs._values[j++];");
      }
      else{
	assert(false);
      }

      return new typeof(return)(num);
    }

    // matrix multiplications
    Matrix!(Row, ColumnR,
	    ReturnTypeOfOperation!(Op, Shape, ShapeR),
	    MatOdr) opBinary(string Op: "*",
			     TypeOfSize ColumnR,
			     MatrixType ShapeR)(in Matrix!(Column, ColumnR, ShapeR, MatOdr) rhs)
    if(ShapeR !is MatrixType.zero){
      enum TypeOfIndex LEN= arrayLength!(Row, ColumnR, ReturnTypeOfOperation!(Op, Shape, ShapeR));
      T[LEN] num= void;

      final switch(ShapeR){
      case MatrixType.dense:
	alias idxMapDest= indexMap!(TemplateArgsOf!(typeof(return)));
	alias idxMapLhs= indexMap!TemplateArgs;
	alias idxMapRhs= indexMap!(TemplateArgsOf!(typeof(rhs)));
	foreach(scope i; 0u..Row)
	  foreach(scope j; 0u..ColumnR){{
	    num[idxMapDest(i, j)]= this._values[idxMapLhs(i, 0u)]*rhs._values[idxMapRhs(j, 0u)];
	    foreach(scope k; 1u..Column) num[idxMapDest(i, j)] += this._values[idxMapLhs(i, k)]*rhs._values[idxMapRhs(k, j)];
	}}
	break;
      case MatrixType.band1:
	alias idxMapLhs= indexMap!TemplateArgs;
	TypeOfIndex k= 0u;
	foreach(scope i; 0u..Row)
	  foreach(scope j; 0u..Column){{
	    num[k++]= this._values[idxMapLhs(i, j)]*rhs._values[j];
	}}
	break;
      case MatrixType.band3:
	size_t idxPrd= 0u, idxLhs= void, idxRhs= void, en= void;
	foreach(scope i; 0u..Row){
	  idxLhs= i*Column;
	  idxRhs= 0u;
	  foreach(scope j; 0u..ColumnR){
	    en= (j == 0u || j == ColumnR-1u)? 2u : 3u;
	    num[idxPrd]= this._values[idxLhs]*rhs._values[idxRhs];
	    foreach(scope k; 1u..en) num[idxPrd] += this._values[idxLhs+k]*rhs._values[idxRhs+2u*k];
	    ++idxPrd;
	    if(j == 0u){
	      idxRhs += 1u;
	    }
	    else{
	      ++idxLhs;
	      idxRhs += 3u;
	    }
	  }
	}
	break;
      case MatrixType.upperTri:
	size_t temp= void;
	alias idxMapDest= indexMap!(TemplateArgsOf!(typeof(return)));
	alias idxMapLhs= indexMap!TemplateArgs;
	foreach(scope i; iota!TypeOfSize(ColumnR)){
	  foreach(scope j; iota!TypeOfSize(Row)){
	    num[idxMapDest(j, i)]= this._values[idxMapLhs(j, 0u)]*rhs._values[i];
	    temp= 0u;
	    foreach(scope k; 1u..i+1u){
	      num[idxMapDest(j, i)] += this._values[idxMapLhs(j, k)]*rhs._values[i-1u+ColumnR+temp];
	      temp += ColumnR-(k+1);
	    }
	  }
	}
	break;
      case MatrixType.lowerTri:
	import sekitk.integers.progression: sumFromZero;

	size_t indexRhs(in size_t i, in size_t k) @safe pure nothrow @nogc const{
	  return i*(i+3u)/2u+sumFromZero(k)-sumFromZero(i);
	}

	size_t idxPrd= void, idxRhs= void, temp= void;
	foreach(scope i; 0u..Column){
	  foreach(scope j; 0u..Row){
	    idxPrd= j*ColumnR+i;
	    idxRhs= indexRhs(i, ColumnR-1u);
	    num[idxPrd]= this._values[Column*(j+1u)-1u]*rhs._values[idxRhs];
	    foreach_reverse(scope k; i..ColumnR-1u) num[idxPrd] += this._values[j*Column+k]*rhs._values[indexRhs(i, k)];
	  }
	}
	break;
      case MatrixType.zero, MatrixType.permutation:
	assert(false);
      }
      return new typeof(return)(num);
    }
    unittest{
      immutable float[21] tempDense1= [2.0, 5, 11, 17, 23, 31, 41,
				       47, 59, 67, 73, 83, 97, 103,
				       109, 127, 137, 149, 157, 167, 173];
      immutable float[35] tempDense2= [3.0, 7, 13, 19, 29,
				       37, 43, 53, 61, 71,
				       79, 89, 101, 107, 113,
				       131, 139, 151, 163, 179,
				       181, 191, 193, 197, 199,
				       211, 223, 227, 229, 233,
				       239, 241, 251, 257, 263];
      immutable float[15] result= [23790.0, 24758, 25736, 26458, 27282,
				   77287, 81283, 85419, 88595, 92397,
				   140369, 148049, 156117, 162397, 169983];

      auto mDense1= new SekiTK!float.Matrix!(3, 7, MatrixType.dense)(tempDense1);
      auto mDense2= new SekiTK!float.Matrix!(7, 5, MatrixType.dense)(tempDense2);
      auto mul= mDense1*mDense2;
      assert(mul._values[] == result[]);
    }
    unittest{
      SekiTK!float.Matrix!(4, 5, MatrixType.dense) lhsDense;
      SekiTK!float.Matrix!(5, 5, MatrixType.band1) rhsDiag;
      SekiTK!float.Matrix!(5, 5, MatrixType.band3) rhsBand3;
      SekiTK!float.Matrix!(5, 5, MatrixType.upperTri) rhsUptri;
      SekiTK!float.Matrix!(5, 5, MatrixType.lowerTri) rhsLotri;
      {
	immutable float[20] tempDense= [2.0, 3, 5, 7, 11,
					13, 17, 19, 23, 29,
					31, 37, 41, 43, 47,
					53, 59, 61, 67, 71];
	lhsDense= new typeof(lhsDense)(tempDense);
	immutable float[5] tempDiag= [73, 79, 83, 89, 97];
	rhsDiag= new typeof(rhsDiag)(tempDiag);
	immutable float[13] tempBand3= [157, 151,
					149, 139, 137,
					131, 127, 113,
					107,101, 89,
					79, 73];
	rhsBand3= new typeof(rhsBand3)(tempBand3);
	immutable float[15] tempUptri=[73.0, 79, 83, 89, 97,
				       101, 103, 107, 109,
				       113, 127, 131,
				       137, 139,
				       149];
	rhsUptri= new typeof(rhsUptri)(tempUptri);
	immutable float[15] tempLotri=[73.0,
				       79, 83,
				       89, 97, 101,
				       103, 107, 109, 113,
				       127, 131, 137, 139, 149];
	rhsLotri= new typeof(rhsLotri)(tempLotri);
      }
      // dense * diagonal
      {
	immutable float[20] result= [146, 237, 415, 623, 1067,
				     949, 1343, 1577, 2047, 2813,
				     2263, 2923, 3403, 3827, 4559,
				     3869, 4661, 5063, 5963, 6887];
	auto mul= lhsDense*rhsDiag;
	assert(mul._values[] == result[]);
      }
      // dense * band3
      {
	immutable float[20] result= [761, 1374, 1795, 2141, 1426,
				     4574, 6815, 7203, 6761, 4164,
				     10380, 15195, 14877, 12689, 7258,
				     17112, 24195, 22999, 19269, 11146];
	auto mul= lhsDense*rhsBand3;
	assert(mul._values[] == result[]);						
      }
      // dense * upper triangular
      {
	immutable float[20] result= [146.0, 461, 1040, 2093, 3788,
				     949, 2744, 4977, 8540, 13121,
				     2263, 6186, 11017, 17816, 25391,
				     3869, 10146, 17369, 27956, 39455];
	auto mul= lhsDense*rhsUptri;
	assert(mul._values[] == result[]);
      }
      // dense * lower triangular
      {
	immutable float[20] result= [2946.0, 2924, 2775, 2320, 1639,
				     10035, 9514, 8399, 6630, 4321,
				     19233, 17806, 15267, 11392, 7003,
				     29877, 27284, 23191, 17440, 10579];
	auto mul= lhsDense*rhsLotri;
	assert(mul._values[] == result[]);
      }
    }	// end of unittest
  }

  const{
    /**
     * ---
     * auto mat= new Matrix!(2, 3, MatrixType.dense, MajorOrder.row);
     * foreach(T elm; mat) ...;
     * ---
     *
     * FIXME
     */
    int opApply(Dg)(scope Dg dg)
    if(ParameterTypeTuple!Dg.length == 1){
      typeof(return) result= 0;

      foreach(scope elm; _values){
	result= dg(elm);
      }
      return result;
    }

    /**
     * ---
     * auto mat= new Matrix!(2, 3, MatrixType.dense, MajorOrder.row);
     * foreach(Vector!2 vec; mat) ...;
     * ---
     *
     * FIXME
     */
    int opApply(Dg)(scope Dg dg)
    if(ParameterTypeTuple!Dg.length == 1){
      typeof(return) result= 0;

      foreach(scope idx; 0..Column){
	result= dg(Vector!Row(sliceColumn(idx)));
      }
      return result;
    }
  }
};


/*************************************************************
 * Rectangular dense matrix
 *************************************************************/
/+
mixin template RectDense(T, real Threshold,
												 TypeOfSize Row,
												 TypeOfSize Column,
												 MajorOrder MatOdr)
if(Row > 0 && Column > 0){
+/
enum string RECT_DENSE= q{
  mixin(DENSE_COMMON);

  version(future){
    void pseudoInverse(){};
  }
};


/*************************************************************
 * Square dense matrix
 *************************************************************/
//mixin template SqDense(T, real Threshold, TypeOfSize Size, MajorOrder MatOdr)
//if(Size > 0){
enum string SQ_DENSE= q{
  mixin(DENSE_COMMON);

private @safe pure:
  /********************************************
   * determinant
   ********************************************/
  T detImpl() nothrow @nogc const{
  import numeric.pseudocmplx: assumeRealNum;
  typeof(return) result= void;

  static if(Size == 1){
    result= _values[0];
  }
  else static if(Size == 2){
    /// Standars: rule of Sarrus
    result= _values[0]*_values[3] -_values[1]*_values[2];
  }
  else static if(Size == 3){
    /// Standards: rule of Sarrus
    result= _values[0]*(_values[4]*_values[8] -_values[5]*_values[7])
      +_values[1]*(_values[5]*_values[6]-_values[3]*_values[8])
      +_values[2]*(_values[3]*_values[7]-_values[4]*_values[6]);
  }
  else{
    /+
     if(!_ev.isNull){
     static if(isComplex!T) result= ev.det;
     else result= assumeRealNum!(CT, Threshold)(ev.det);
     }
     else if(!resultLU.isNull) result= resultLU.get.det;	// LU decomposition is obtained
     else{
     auto temp= ForwardElimination!(Size, Shape, MatOdr).transform(this.opCast!(T[Column][Row]));	// forward elimination
     result= temp.det;
     }+/
    result= T.nan;
  }
  return result;
}

  /********************************************
   * inverse
   ********************************************/
  static if(Size == 1){
    TypeOfInternalArray inverseImpl() @safe pure nothrow @nogc const{
      return [Identity!(T, "*")/_values[0]];
    }
  }
  else static if(Size == 2){
    TypeOfInternalArray inverseImpl() @safe pure nothrow @nogc const{
      typeof(return) result= [_values[3], -_values[1], -_values[2], _values[0]];
      result[] /= this.detImpl;
      return result;
    }
  }
  else{
    TypeOfInternalArray inverseImpl() @safe pure{
      typeof(return) num= void;
/+
  if(resultLU.isNull) this.decompose!(DecompScheme.luWithPartialPivDoolittle);
    alias TypCurr= resultLU.get.type;
    if(resultLU.get.get!TypCurr.isInvertible){
    num= resultLU.get.get!TypCurr.inverseArray;
  }+/
      return num;
    }
  }

	/+
a11, a12
a21, a22

1,   0| u11 u12| u11      u12
l21, 1|   0 u22| l21*u11  l21*u12+u22

u11= a11
u12= a12
l21*u11= a21
    l21= a21/a11
l21*u12+u22= a22
u22= a22-a21*a12/a11
+/
};
