#include "vec3.h"

Vec3::Vec3() {
	p = new double[3];
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;
}

Vec3::Vec3(double x_, double y_, double z_) {
	p = new double[3];
	p[0] = x_;
	p[1] = y_;
	p[2] = z_;
}

Vec3::Vec3(double p_[3]) {
	p = new double[3];
	p[0] = p_[0];
	p[1] = p_[1];
	p[2] = p_[2];
}

Vec3::~Vec3() {
	delete[] p;
}

double * Vec3::data() const {
	return p;
}

Vec3 Vec3::operator+(const Vec3 & v) const {
	return Vec3(p[0] + v.p[0], p[1] + v.p[1], p[2] + v.p[2]);
}

Vec3 Vec3::operator-(const Vec3 & v) const {
	return Vec3(p[0] - v.p[0], p[1] - v.p[1], p[2] - v.p[2]);
}

Vec3 Vec3::operator*(double c) const {
	return Vec3(p[0] * c, p[1] * c, p[2] * c);
}

double Vec3::operator*(const Vec3 & v) const {
	return p[0] * v.p[0] + p[1] * v.p[1] + p[2] * v.p[2];
}

Vec3 Vec3::operator/(double c) const {
	double c_ = 1.0 / c;
	return Vec3(p[0] * c_, p[1] * c_, p[2] * c_);
}
