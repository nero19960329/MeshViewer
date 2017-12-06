#pragma once

class Vec3 {
private:
	double * p;

public:
	Vec3();
	Vec3(double x_, double y_, double z_);
	Vec3(double p_[3]);
	~Vec3();

	double * data() const;

	Vec3 operator + (const Vec3 & v) const;
	Vec3 operator - (const Vec3 & v) const;
	Vec3 operator * (double c) const;
	double operator * (const Vec3 & v) const;
	Vec3 operator / (double c) const;
};