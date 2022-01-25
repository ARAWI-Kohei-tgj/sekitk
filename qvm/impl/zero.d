module sekitk.qvm.impl.zero;
/+
import std.traits: isFloatingPoint;
import numeric.pseudocmplx: isComplex;
mixin template MatrixImplZero(T, T Threshold)
if(isFloatingPoint!T || isComplex!T){
	import sekitk.base: TypeOfSize, MatrixType, MajorOrder, matrixConstraint;
+/
enum string MATRIX_IMPL_ZERO= q{
  /************************************************************
   * zero matrix
   ************************************************************/
  class Matrix(TypeOfSize Row, TypeOfSize Column,
	       MatrixType Shape: MatrixType.zero,
	       MajorOrder MatOdr= MajorOrder.row)
  if(matrixConstraint!(Row, Column, Shape)){
    import std.traits: Unqual;

    private{
      alias TypeOfThis= typeof(this);
    }

    @safe pure nothrow @nogc{
      this(){}
      this(in TypeOfThis other){}
    }

    @safe pure nothrow const{
      // sign
      TypeOfThis opUnary(string Op)() @nogc
      if(isPlusOrMinusSign!Op){
	return new typeof(return)();
      }

      // add & sub
      Matrix!(Row, Column,
	      ShapeR,
	      MatOdr) opBinary(string Op,
			       MatrixType ShapeR)(in Matrix!(Row, Column,
							     ShapeR,
							     MatOdr) rhs) @nogc
      if(isPlusOrMinusSign!Op){
	return new typeof(return)(rhs);
      }

      /****************************************
       * O * scalar, O / scalar
       ****************************************/
      TypeOfThis opBinary(string Op, TypeR)(in TypeR rhs) @nogc
      if((Op == "*" || Op == "/") && is(TypeR: T)){
	return new typeof(return)();
      }

      // scalar * matrix
      // O(m*n)
      TypeOfThis opBinaryRight(string OP: "*")(in T lhs){
	return new typeof(return)();
      }

      // matrix * vector
      Matrix!(Row, 1, Shape, MatOdr) opBinary(string OP: "*")(in Vector!Column rhs){
	return new typeof(return)();
      }

      /****************************************
       * O * Matrix
       ****************************************/
      Matrix!(Row, ColumnR,
	      Shape, MatOdr) opBinary(string OP: "*",
				      TypeOfSize ColumnR,
				      MatrixType ShapeR)(in Matrix!(Column, ColumnR,
								    ShapeR, MatOdr) rhs){
	return new typeof(return)();
      }
    }
  }
};
