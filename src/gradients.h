#include "largeVis.h"

class Gradient {
protected:
	const distancetype gamma;
	distancetype cap;
	const dimidxtype D;
	Gradient(const distancetype& g, const dimidxtype& d);
	virtual void _positiveGradient(const distancetype& dist_squared,
                                 coordinatetype* holder) const = 0;
	virtual void _negativeGradient(const distancetype& dist_squared,
                                 coordinatetype* holder) const = 0;
	void multModify(coordinatetype *col, const coordinatetype& adj) const;
	void multModifyPos(coordinatetype *col, const coordinatetype& adj) const;
	coordinatetype clamp(const coordinatetype& val) const;
	distancetype distAndVector(const coordinatetype *x_i,
                             const coordinatetype *x_j,
                             coordinatetype *output) const;

public:
	virtual ~Gradient();
	void positiveGradient(const coordinatetype* i,
	                      const coordinatetype* j,
	                      coordinatetype* holder) const;
	void negativeGradient(const coordinatetype* i,
                        const coordinatetype* k,
                        coordinatetype* holder) const;
};

class AlphaGradient: public Gradient {
	const coordinatetype alpha;
	const coordinatetype twoalpha;
protected:
	const coordinatetype alphagamma;
	virtual void _positiveGradient(const double& dist_squared,
                                coordinatetype* holder) const;
	virtual void _negativeGradient(const double& dist_squared,
                                coordinatetype* holder) const;
public:
	AlphaGradient(const distancetype& a,
                const distancetype& g,
                const dimidxtype& D);
};


class AlphaOneGradient: public AlphaGradient {
public:
	AlphaOneGradient(const distancetype& g,
                   const dimidxtype& d);
protected:
	virtual void _positiveGradient(const distancetype& dist_squared,
                                coordinatetype* holder) const;
	virtual void _negativeGradient(const distancetype& dist_squared,
                                coordinatetype* holder) const;
};

class ExpGradient: public Gradient {
public:
	const coordinatetype gammagamma;
	ExpGradient(const distancetype& g, const dimidxtype& d);
protected:
	virtual void _positiveGradient(const distancetype& dist_squared,
                                coordinatetype* holder) const;
	virtual void _negativeGradient(const distancetype& dist_squared,
                                coordinatetype* holder) const;
};
