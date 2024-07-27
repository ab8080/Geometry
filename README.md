# Geometry

In this task, you need to implement a set of classes for solving geometric problems on a plane. All coordinates are assumed to be integers.

## Vector Class

For a vector on the plane, you need to implement the following methods and overload the corresponding operators:
- Default constructor (creates a zero vector)
- Constructor from two integer variables (creates a vector with corresponding coordinates)
- Operator `*` for dot product
- Operator for cross product
- Operators for addition/subtraction with another vector, accordingly implement operators `+=` and `-=`
- Operator for multiplication by a number (make it so that both vector * number and number * vector are allowed), implement operator `*=`
- Unary minus operator to get a vector in the opposite direction

## Shape Classes

Create a set of classes for shapes that inherit from the abstract class `IShape` to work with 2D primitives:
- Point
- Segment
- Line
- Ray
- Polygon
- Circle

In the base class, provide the following methods:
- `void Move(const Vector&)` — moves the shape by the corresponding vector
- `bool ContainsPoint(const Point&)` — checks if the shape contains the point
- `bool CrossSegment(const Segment&)` — checks if the segment intersects the shape
- `IShape* Clone()` — returns a pointer to a copy of the shape
- `std::string ToString()` — string representation of the shape (format in examples)

In derived classes, you need to implement these methods.
