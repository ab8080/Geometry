#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace Geometry {

class Vector;
class IShape;
class Point;
class Segment;
class Ray;
class Line;
class Circle;
class Polygon;

class Vector {
 public:
  Vector();
  Vector(const double& x1, const double& y1);
  Vector(const Point& p1, const Point& p2);

  [[nodiscard]] double Len() const;

  double operator*(Vector& b) const;
  double operator^(Vector& b) const;
  Vector& operator+=(const Vector& b);
  Vector& operator-=(const Vector& b);
  Vector& operator*=(const double& n);
  Vector operator+(const Vector& b) const;
  Vector operator-(const Vector& b) const;

  virtual void Print() const;
  double TriangleArea(Vector& b) const;
  [[nodiscard]] Point GetXY() const;

 private:
  double x_ = 0;
  double y_ = 0;
};

class IShape : public Vector {
 public:
  virtual IShape& Move(const Point&) = 0;
  [[nodiscard]] virtual bool ContainsPoint(const Point&) const = 0;
  [[nodiscard]] virtual bool CrossesSegment(const Segment&) const = 0;
  [[nodiscard]] virtual IShape* Clone() const = 0;
  [[nodiscard]] virtual std::string ToString() const = 0;
  virtual ~IShape() = default;
};

class Point : public IShape {
 public:
  Point() = default;
  Point(const double& a, const double& b);
  Point(const int& a, const int& b);

  Point& operator+=(const Point& b);
  Point& operator-=(const Point& b);
  Point operator-(const Point& b) const;
  Point& operator=(const Point& b);

  [[nodiscard]] std::string ToString() const override;
  Point& Move(const Point& v) override;
  [[nodiscard]] bool ContainsPoint(const Point& p) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& p) const override;
  [[nodiscard]] Point* Clone() const override;
  ~Point() override = default;
  [[nodiscard]] std::pair<double, double> Getxy() const;

 private:
  double x1_ = 0;
  double y1_ = 0;
  double x2_ = 0;
  double y2_ = 0;
};

class Line : public IShape {
 public:
  Line() = default;
  Line(const double& a, const double& b, const double& c);
  Line(const Point& p1, const Point& p2);
  [[nodiscard]] bool Crossing(const Line& b) const;
  [[nodiscard]] Point CrossPoint(const Line& l) const;
  [[nodiscard]] double Distance(const Line& l) const;
  [[nodiscard]] double Distance(const Point& p) const;
  void Print() const override;
  [[nodiscard]] Line Copy() const;
  Line& Move(const Point& v) override;
  [[nodiscard]] bool ContainsPoint(const Point& p) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& s) const override;
  [[nodiscard]] Line* Clone() const override;
  [[nodiscard]] std::string ToString() const override;
  ~Line() override = default;

 private:
  double x_ = 0;
  double y_ = 0;
  double a_ = 0;
  double b_ = 0;
  double c_ = 0;
  Point point_a_ = {0, 0};
  Point point_b_ = {0, 0};
};

class Ray : public IShape {
 public:
  Ray() = default;
  Ray(const Point& p, const Point& l);
  [[nodiscard]] bool ContainsPoint(const Point& p) const override;
  [[nodiscard]] double Distance(const Point& p) const;
  [[nodiscard]] bool CrossesSegment(const Segment& s) const override;
  Ray& Move(const Point& v) override;
  [[nodiscard]] Ray* Clone() const override;
  [[nodiscard]] std::string ToString() const override;
  ~Ray() override = default;

 private:
  Point start_;
  Point other_;
};

class Segment : public IShape {
 public:
  Segment() = default;
  Segment(const Point& p, const Point& l);
  [[nodiscard]] bool ContainsPoint(const Point& p) const override;
  [[nodiscard]] bool Crossing(const Segment& s) const;
  [[nodiscard]] bool CrossesSegment(const Segment& s) const override;
  [[nodiscard]] double Distance(const Point& p) const;
  Segment& Move(const Point& v) override;
  [[nodiscard]] Segment* Clone() const override;
  [[nodiscard]] std::string ToString() const override;
  [[nodiscard]] std::pair<Point, Point> GetStartEnd() const;
  ~Segment() override = default;

 private:
  Point start_;
  Point end_;
};

class Polygon : public IShape {
 public:
  Polygon() = default;
  explicit Polygon(const std::vector<Point>& arr);
  bool Convex();
  [[nodiscard]] bool ContainsPoint(const Point& p) const override;
  Polygon& Move(const Point& v) override;
  [[nodiscard]] bool CrossesSegment(const Segment& s) const override;
  [[nodiscard]] Polygon* Clone() const override;
  [[nodiscard]] std::string ToString() const override;
  ~Polygon() override = default;

 private:
  std::vector<Point> points_;
};

class Circle : public IShape {
 public:
  Circle() = default;
  Circle(const Point& centre, int radius);
  [[nodiscard]] bool ContainsPoint(const Point& p) const override;
  [[nodiscard]] bool ContainsPointInside(const Point& p) const;
  Circle& Move(const Point& v) override;
  [[nodiscard]] bool CrossesSegment(const Segment& s) const override;
  [[nodiscard]] Circle* Clone() const override;
  [[nodiscard]] std::string ToString() const override;
  ~Circle() override = default;

 private:
  Point centre_;
  int radius_ = 0;
};

Point::Point(const int& a, const int& b) {
  x1_ = a;
  y1_ = b;
}

Point::Point(const double& a, const double& b) {
  x1_ = a;
  y1_ = b;
}

Point Point::operator-(const Point& b) const {
  Point copy(x1_ - b.x1_, y1_ - b.y1_);
  return copy;
}

Point& Point::operator=(const Point& b) {
  x1_ = b.x1_;
  y1_ = b.y1_;
  return *this;
}

Point& Point::operator+=(const Point& b) {
  x1_ += b.x1_;
  y1_ += b.y1_;
  return *this;
}

Point& Point::operator-=(const Point& b) {
  x1_ -= b.x1_;
  y1_ -= b.y1_;
  return *this;
}

std::string Point::ToString() const {
  std::stringstream s;
  s << "Point(" << x1_ << ", " << y1_ << ")";
  return s.str();
}

Point& Point::Move(const Point& v) {
  *this += v;
  return *this;
}

bool Point::ContainsPoint(const Point& p) const {
  return (x1_ == p.x1_ && y1_ == p.y1_);
}

bool Point::CrossesSegment(const Segment& p) const {
  return p.ContainsPoint(*this);
}

Point* Point::Clone() const {
  auto* new_point = new Point(x1_, y1_);
  return new_point;
}

std::pair<double, double> Point::Getxy() const {
  return std::pair<double, double>(x1_, y1_);
}

Circle::Circle(const Point& centre, int radius) {
  centre_ = centre;
  radius_ = radius;
}

bool Circle::ContainsPoint(const Point& p) const {
  Vector ab(centre_, p);
  if (ab.Len() == radius_) {
    return true;
  }
  return ab.Len() < radius_;
}

Circle& Circle::Move(const Point& v) {
  centre_ += v;
  return *this;
}

bool Circle::CrossesSegment(const Segment& s) const {
  double d1 = s.Distance(centre_);
  auto points = s.GetStartEnd();
  if (d1 < radius_) {
    return !(this->ContainsPointInside(points.first) &&
             this->ContainsPointInside(points.second));
  }
  return d1 <= radius_;
}

Circle* Circle::Clone() const {
  auto* new_circle = new Circle;
  *new_circle = Circle(centre_, radius_);
  return new_circle;
}

std::string Circle::ToString() const {
  std::stringstream s;
  s << "Circle("
    << "Point(" << centre_.Getxy().first << ", " << centre_.Getxy().second
    << ")"
    << ", " << radius_ << ")";
  return s.str();
}

bool Circle::ContainsPointInside(const Point& p) const {
  Vector ab(centre_, p);
  if (ab.Len() == radius_) {
    return false;
  }
  return ab.Len() < radius_;
}

Polygon::Polygon(const std::vector<Point>& arr) {
  for (const auto& i : arr) {
    points_.push_back(i);
  }
}

bool Polygon::Convex() {
  uint64_t size = points_.size();
  std::vector<Vector> vectors(size);
  vectors[size - 1] = Vector(points_[size - 1], points_[0]);
  for (uint64_t i = 0; i < size - 1; ++i) {
    vectors[i] = Vector(points_[i], points_[i + 1]);
  }
  double sign = (vectors[0] ^ vectors[1]);
  for (uint64_t i = 2; i < size; ++i) {
    if ((vectors[i - 1] ^ vectors[i]) * sign < 0) {
      return false;
    }
  }
  return (vectors[size - 1] ^ vectors[0]) * sign >= 0;
}

Polygon& Polygon::Move(const Point& v) {
  for (auto& point : points_) {
    point += v;
  }
  return *this;
}

bool Polygon::CrossesSegment(const Segment& s) const {
  for (uint64_t i = 0; i < points_.size(); ++i) {
    uint64_t last = (i + 1) % points_.size();
    Line l(points_[i], points_[last]);
    if (l.CrossesSegment(s)) {
      return true;
    }
  }
  return false;
}

Polygon* Polygon::Clone() const {
  auto* new_poly = new Polygon;
  *new_poly = Polygon(points_);
  return new_poly;
}

std::string Polygon::ToString() const {
  std::stringstream s;
  s << "Polygon(";
  for (uint64_t i = 0; i < points_.size() - 1; ++i) {
    std::string p = points_[i].ToString();
    s << p << ", ";
  }
  s << points_[points_.size() - 1].ToString();
  s << ")";
  return s.str();
}

bool Polygon::ContainsPoint(const Point& p) const {
  for (uint64_t i = 0; i < points_.size(); ++i) {
    uint64_t last = (i + 1) % points_.size();
    Line l(points_[i], points_[last]);
    if (l.ContainsPoint(p)) {
      return true;
    }
  }
  // we've checked if borders contain point
  const double kBigNum = 100001;
  Point far(kBigNum, p.Getxy().second + 1);

  int counter = 0;
  for (uint64_t i = 0; i < points_.size(); ++i) {
    uint64_t last = (i + 1) % points_.size();
    Segment s1(p, far);
    Segment s2(points_[i], points_[last]);
    if (s1.Crossing(s2)) {
      ++counter;
    }
  }
  return ((counter % 2) == 1);
}

Segment::Segment(const Point& p, const Point& l) {
  start_ = p;
  end_ = l;
}

bool Segment::ContainsPoint(const Point& p) const {
  const Point kZero(0, 0);
  Vector ap = Vector(kZero, p - start_);
  Vector curr(start_, end_);
  if ((curr ^ ap) == 0) {
    double x = p.Getxy().first;
    double y = p.Getxy().second;
    double start_x = start_.Getxy().first;
    double start_y = start_.Getxy().second;
    double end_x = end_.Getxy().first;
    double end_y = end_.Getxy().second;
    return (x <= std::max(start_x, end_x) && x >= std::min(start_x, end_x) &&
            y <= std::max(start_y, end_y) && y >= std::min(start_y, end_y));
  }
  return false;
}

bool Segment::Crossing(const Segment& s) const {
  Vector ab = Vector(start_, end_);
  Vector ac = Vector(start_, s.start_);
  Vector ad = Vector(start_, s.end_);
  Vector cd = Vector(s.start_, s.end_);
  Vector cb = Vector(s.start_, end_);
  Vector db = Vector(s.end_, end_);
  return (((ab ^ ac) * (ab ^ ad) < 0 && (ac ^ cd) * (cb ^ cd) > 0) ||
          ((ab ^ ac) == 0 && (ac * cb) >= 0) ||
          ((ab ^ ad) == 0 && (ad * db) >= 0) ||
          ((ac ^ cd) == 0 && ac * ad <= 0) || ((cb ^ cd) == 0 && cb * db <= 0));
  // k
}

double Segment::Distance(const Point& p) const {
  Segment s(start_, end_);
  if (s.ContainsPoint(p)) {
    return 0;
  }
  Vector ab(start_, end_);
  Vector ap(start_, p);
  Vector bp(end_, p);
  if ((ab * ap) > 0 && (bp * ab) <= 0) {
    Line l(start_, end_);
    return l.Distance(p);
  }
  return std::min(ap.Len(), bp.Len());
}

bool Segment::CrossesSegment(const Segment& s) const {
  return this->Crossing(s);
}

Segment& Segment::Move(const Point& v) {
  start_ += v;
  end_ += v;
  return *this;
}

Segment* Segment::Clone() const {
  auto* new_seg = new Segment;
  *new_seg = Segment(start_, end_);
  return new_seg;
}

std::string Segment::ToString() const {
  std::stringstream s;
  std::string p1 = start_.ToString();
  std::string p2 = end_.ToString();
  s << "Segment(" << p1 << ", " << p2 << ")";
  return s.str();
}

std::pair<Point, Point> Segment::GetStartEnd() const { return {start_, end_}; }

Ray::Ray(const Point& p, const Point& l) {
  start_ = p;
  other_ = l;
}

bool Ray::ContainsPoint(const Point& p) const {
  Vector v1(start_, other_);
  Vector v2(start_, p);
  Line l(start_, other_);
  return (v1 * v2) >= 0 && l.ContainsPoint(p);
}

double Ray::Distance(const Point& p) const {
  Ray r(start_, other_);
  if (r.ContainsPoint(p)) {
    return 0;
  }
  Vector ab(start_, other_);
  Vector ap(start_, p);
  if ((ab * ap) > 0) {
    Line l(start_, other_);
    return l.Distance(p);
  }
  return ap.Len();
}

bool Ray::CrossesSegment(const Segment& s) const {
  Vector v(start_, other_);
  Vector v1(start_, s.GetStartEnd().first);
  Vector v2(start_, s.GetStartEnd().second);
  if ((v1 ^ v2) == 0) {
    if ((v1 ^ v) == 0) {
      return ((v1 * v) > 0 || (v2 * v) > 0 || (v1 * v) == 0 || (v2 * v) == 0);
    }
    return (v1 * v2) <= 0;
  }
  if (((v1 ^ v) > 0 && (v ^ v2) > 0 && (v1 ^ v2) > 0) ||
      ((v1 ^ v) < 0 && (v ^ v2) < 0 && (v1 ^ v2) < 0)) {
    return true;
  }
  return (((v ^ v2) == 0 && (v * v2) > 0) || ((v1 ^ v) == 0 && (v1 * v) > 0));
}

Ray& Ray::Move(const Point& v) {
  start_ += v;
  other_ += v;
  return *this;
}

Ray* Ray::Clone() const {
  auto* new_ray = new Ray;
  *new_ray = Ray(start_, other_);
  return new_ray;
}

std::string Ray::ToString() const {
  double x = other_.Getxy().first - start_.Getxy().first;
  double y = other_.Getxy().second - start_.Getxy().second;
  std::stringstream s;
  s << "Ray("
    << "Point(" << start_.Getxy().first << ", " << start_.Getxy().second << ")"
    << ", Vector(" << x << ", " << y << "))";
  return s.str();
}

Line::Line(const double& a, const double& b, const double& c) {
  a_ = a;
  b_ = b;
  c_ = c;
  if (a == 0) {
    x_ = 1;
    y_ = 0;
  } else if (b == 0) {
    x_ = 0;
    y_ = 1;
  } else {
    x_ = 1;
    y_ = -a / b;
  }
}

bool Line::Crossing(const Line& b) const {
  if (this->Distance(b) == 0) {
    return true;
  }
  Vector v1(x_, y_);
  Vector v2(b.x_, b.y_);
  return (v1 ^ v2) != 0;
}

void Line::Print() const { std::cout << x_ << ' ' << y_ << '\n'; }

double Line::Distance(const Line& l) const {
  double d;
  if (l.a_ == 0 && a_ != 0) {
    d = (l.c_ - (l.a_ / a_) * c_) / sqrt(l.a_ * l.a_ + l.b_ * l.b_);
  } else if (l.a_ == 0 && a_ == 0) {
    d = -c_ / b_ - l.c_ / l.b_;
  } else {
    d = (c_ - (a_ / l.a_) * l.c_) / sqrt(a_ * a_ + b_ * b_);
  }
  if (d < 0) {
    return -d;
  }
  return d;
}

double Line::Distance(const Point& p) const {
  if (this->ContainsPoint(p)) {
    return 0;
  }
  double d = (a_ * p.Getxy().first + b_ * p.Getxy().second + c_) /
             sqrt(a_ * a_ + b_ * b_);
  if (d < 0) {
    return -d;
  }
  return d;
}

Point Line::CrossPoint(const Line& l) const {
  if (a_ != 0) {
    double part1 = c_ * l.a_ / a_ - l.c_;
    double part2 = -b_ * l.a_ / a_ + l.b_;
    double y = part1 / part2;
    double x = -(b_ * y + c_) / a_;
    if (x == -0) {
      x = 0;
    }
    if (y == -0) {
      y = 0;
    }
    Point p = {x, y};
    return p;
  }
  double y = -c_ / b_;
  double x = -(l.b_ * y + l.c_) / l.a_;
  Point p = {x, y};
  return p;
}

Line::Line(const Point& p1, const Point& p2) {
  a_ = p2.Getxy().second - p1.Getxy().second;
  b_ = -p2.Getxy().first + p1.Getxy().first;
  c_ = p2.Getxy().first * p1.Getxy().second -
       p1.Getxy().first * p2.Getxy().second;
  x_ = p2.Getxy().first - p1.Getxy().first;
  y_ = p2.Getxy().second - p1.Getxy().second;
  point_a_ = p1;
  point_b_ = p2;
}

bool Line::ContainsPoint(const Point& p) const {
  return (a_ * p.Getxy().first + b_ * p.Getxy().second + c_ == 0);
}

Line Line::Copy() const {
  Line copy = *this;
  return copy;
}

bool Line::CrossesSegment(const Segment& s) const {
  Ray r1(point_a_, point_b_);
  Ray r2(point_b_, point_a_);
  return (r1.CrossesSegment(s) || r2.CrossesSegment(s));
}

Line& Line::Move(const Point& v) {
  c_ = c_ - a_ * v.Getxy().first - b_ * v.Getxy().second;
  return *this;
}

Line* Line::Clone() const {
  auto* new_line = new Line;
  *new_line = Line(a_, b_, c_);
  return new_line;
}

std::string Line::ToString() const {
  std::stringstream s;
  if (c_ == -0) {
    s << "Line(" << a_ << ", " << b_ << ", " << 0 << ")";
  } else {
    s << "Line(" << a_ << ", " << b_ << ", " << c_ << ")";
  }
  return s.str();
}

Vector::Vector() {
  x_ = 0;
  y_ = 0;
}

Vector::Vector(const double& x1, const double& y1) {
  x_ = x1;
  y_ = y1;
}

Vector::Vector(const Point& p1, const Point& p2) {
  x_ = p2.Getxy().first - p1.Getxy().first;
  y_ = p2.Getxy().second - p1.Getxy().second;
}

double Vector::operator*(Vector& b) const { return (x_ * b.x_ + b.y_ * y_); }
double Vector::operator^(Vector& b) const { return (x_ * b.y_ - b.x_ * y_); }
Vector& Vector::operator+=(const Vector& b) {
  this->x_ += b.x_;
  this->y_ += b.y_;
  return *this;
}

Vector Vector::operator+(const Vector& b) const {
  Vector copy = Vector(x_, y_);
  copy += b;
  return copy;
}

Vector& Vector::operator-=(const Vector& b) {
  this->x_ -= b.x_;
  this->y_ -= b.y_;
  return *this;
}

Vector Vector::operator-(const Vector& b) const {
  Vector copy = Vector(x_, y_);
  copy -= b;
  return copy;
}

double Vector::Len() const { return sqrt(x_ * x_ + y_ * y_); }

void Vector::Print() const {
  const int kPrecision = 9;
  std::cout.precision(kPrecision);
  std::cout << x_ << ' ' << y_ << '\n';
}

double Vector::TriangleArea(Vector& b) const {
  double area = (*this ^ b);
  if (area < 0) {
    area *= -1;
  }
  return area / 2;
}

Vector& Vector::operator*=(const double& n) {
  x_ *= n;
  y_ *= n;
  return *this;
}

Point Vector::GetXY() const {
  Point p(x_, y_);
  return p;
}
}  // namespace Geometry

template <class SmartPtrT>
void Delete(const SmartPtrT& unused) {}

template <class T>
void Delete(T* ptr) {
  delete ptr;
}

void CheckFunctions(const Geometry::IShape* shape_ptr,
                    const Geometry::Point& point_a,
                    const Geometry::Point& point_b) {
  std::cout << "Given shape "
            << (shape_ptr->ContainsPoint(point_a) ? "contains"
                                                  : "does not contain")
            << " point A\n";

  const auto kKSegmentAb = Geometry::Segment(point_a, point_b);
  std::cout << "Given shape "
            << (shape_ptr->CrossesSegment(kKSegmentAb) ? "crosses"
                                                       : "does not cross")
            << " segment AB\n";

  const auto kKVectorAb = point_b - point_a;
  auto* const kKClonedShapePtr =
      shape_ptr->Clone();  // may return either raw or smart pointer
  std::cout << kKClonedShapePtr->Move(kKVectorAb).ToString();

  Delete(kKClonedShapePtr);  // raw pointer compatibility
}
int main() {
  std::unique_ptr<Geometry::IShape> shape_ptr;

  std::string command;
  std::cin >> command;

  int x = 0;
  int y = 0;
  int a = 0;
  int b = 0;

  if (command == "point") {
    std::cin >> x >> y;
    shape_ptr = std::make_unique<Geometry::Point>(x, y);
  } else if (command == "segment") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Segment>(Geometry::Point(x, y),
                                                    Geometry::Point(a, b));
  } else if (command == "ray") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Ray>(Geometry::Point(x, y),
                                                Geometry::Point(a, b));
  } else if (command == "line") {
    std::cin >> x >> y >> a >> b;
    shape_ptr = std::make_unique<Geometry::Line>(Geometry::Point(x, y),
                                                 Geometry::Point(a, b));
  } else if (command == "polygon") {
    size_t n_points = 0;
    std::cin >> n_points;
    std::vector<Geometry::Point> points;
    points.reserve(n_points);
    for (size_t i = 0; i < n_points; ++i) {
      std::cin >> x >> y;
      points.emplace_back(x, y);
    }
    shape_ptr = std::make_unique<Geometry::Polygon>(std::move(points));
  } else if (command == "circle") {
    std::cin >> x >> y;
    const auto kCenter = Geometry::Point(x, y);
    uint64_t radius = 0;
    std::cin >> radius;
    shape_ptr = std::make_unique<Geometry::Circle>(kCenter, radius);
  } else {
    std::cerr << "Undefined command" << std::endl;
    return 1;
  }

  std::cin >> x >> y;
  const auto kKPointA = Geometry::Point(x, y);
  std::cin >> x >> y;
  const auto kKPointB = Geometry::Point(x, y);

  CheckFunctions(shape_ptr.get(), kKPointA, kKPointB);
  return 0;
}
