// Copyright 2014 Sebastian Boettcher


class Point2D {
 public:
  Point2D();
  Point2D(float x, float y);
  ~Point2D();

  void setX(const float& x);
  void setY(const float& y);
  void setXY(const float& x, const float& y);
  float getX() const;
  float getY() const;
 private:
  float m_x;
  float m_y;
};

Point2D::Point2D()
  : m_x(0.0), m_y(0.0) {
}

Point2D::Point2D(float x, float y)
  : m_x(0.0), m_y(0.0) {
  m_x = x;
  m_y = y;
}

Point2D::~Point2D() {
}



void Point2D::setX(const float& x) {
  m_x = x;
}
void Point2D::setY(const float& y) {
  m_y = y;
}
void Point2D::setXY(const float& x, const float& y) {
  m_x = x;
  m_y = y;
}

float Point2D::getX() const {
  return m_x;
}

float Point2D::getY() const {
  return m_y;
}


