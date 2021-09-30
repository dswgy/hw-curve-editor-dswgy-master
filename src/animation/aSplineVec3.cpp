#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen\Dense>

#pragma warning(disable:4018)
#pragma warning(disable:4244)


ASplineVec3::ASplineVec3() : mInterpolator(new ABernsteinInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
    if (mInterpolator) delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
    mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
    return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
    mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
    return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
    double fps = getFramerate();

	if (mInterpolator) { delete mInterpolator; }
    switch (type)
    {
	case LINEAR: mInterpolator = new ALinearInterpolatorVec3(); break;
	case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3(); break;
	case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3(); break;
	case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3(); break;
	case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3(); break;
	case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3(); break;
	case LINEAR_EULER: mInterpolator = new AEulerLinearInterpolatorVec3(); break;
	case CUBIC_EULER: mInterpolator = new AEulerCubicInterpolatorVec3(); break;
    };
    
    mInterpolator->setFramerate(fps);
    computeControlPoints();
    cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
    return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys[keyID].second = value;
    computeControlPoints();
    cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
    assert(ID >= 0 && ID < mCtrlPoints.size()+2);
    if (ID == 0)
    {
        mStartPoint = value;
        computeControlPoints(false);
    }
    else if (ID == mCtrlPoints.size() + 1)
    {
        mEndPoint = value;
		computeControlPoints(false);
    }
    else mCtrlPoints[ID-1] = value;
    cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
    mKeys.push_back(Key(time, value));

    if (updateCurve)
    {
        computeControlPoints();
        cacheCurve();
    }
}

int ASplineVec3::insertKey(double time, const vec3& value, bool updateCurve)
{
	if (mKeys.size() == 0)
	{
		appendKey(time, value, updateCurve);
		return 0;
	}

	for (int i = 0; i < mKeys.size(); ++i)
	{
		assert(time != mKeys[i].first);
		if (time < mKeys[i].first)
		{
			mKeys.insert(mKeys.begin() + i, Key(time, value));
			if (updateCurve)
			{
				computeControlPoints();
				cacheCurve();
			}
			return i;
		}
	}

	// Append at the end of the curve
	appendKey(time, value, updateCurve);
	return mKeys.size() - 1;
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
    if (mKeys.size() == 0)
    {
        appendKey(0, value, updateCurve);
    }
    else
    {
        double lastT = mKeys[mKeys.size() - 1].first;
        appendKey(lastT + 1, value, updateCurve);
    }
}

void ASplineVec3::deleteKey(int keyID)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys.erase(mKeys.begin() + keyID);
    computeControlPoints();
    cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID) const
{
    assert(keyID >= 0 && keyID < mKeys.size());
    return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
    return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID) const
{
    assert(ID >= 0 && ID < mCtrlPoints.size()+2);
    if (ID == 0) return mStartPoint;
    else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
    else return mCtrlPoints[ID-1];
}

int ASplineVec3::getNumControlPoints() const
{
    return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
    mKeys.clear();
}

double ASplineVec3::getDuration() const 
{
    return mKeys.size() == 0 ? 0 : mKeys[mKeys.size()-1].first;
}

double ASplineVec3::getNormalizedTime(double t) const 
{
    return (t / getDuration());
}

double ASplineVec3::getKeyTime(int keyID) const
{
	assert(keyID >= 0 && keyID < mKeys.size());
	return mKeys[keyID].first;
}

vec3 ASplineVec3::getValue(double t) const
{
    if (mCachedCurve.size() == 0 || mKeys.size() == 0) return vec3();
	if (t < mKeys[0].first)
		return mCachedCurve[0];
	else
		t -= mKeys[0].first;

    double dt = mInterpolator->getDeltaTime();
    int rawi = (int)(t / dt); // assumes uniform spacing
    double frac = (t - rawi*dt) / dt;

	int i = mLooping? rawi % mCachedCurve.size() : std::min<int>(rawi, mCachedCurve.size() - 1);
	int inext = mLooping ? (i + 1) % mCachedCurve.size() : std::min<int>(i + 1, mCachedCurve.size() - 1);

    vec3 v1 = mCachedCurve[i];
    vec3 v2 = mCachedCurve[inext];
    vec3 v = v1*(1 - frac) + v2 * frac;
    return v;
}

void ASplineVec3::cacheCurve()
{
    mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints(bool updateEndPoints)
{
	if (mKeys.size() >= 2 && updateEndPoints)
	{
		int totalPoints = mKeys.size();

		//If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
		//They lie on the tangent of the first and last interpolation points.
		vec3 tmp = mKeys[0].second - mKeys[1].second;
		double n = tmp.Length();
		mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25; // distance to endpoint is 25% of distance between first 2 points

		tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
		n = tmp.Length();
		mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
	}
	//ifNatural = false;
	mInterpolator->ifNatural = ifNatural;
    mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

vec3* ASplineVec3::getCachedCurveData()
{
	return mCachedCurve.data();
}

vec3 * ASplineVec3::getControlPointsData()
{
	return mCtrlPoints.data();
}

int ASplineVec3::getNumCurveSegments() const
{
    return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
    return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
    mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
    return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
    return mDt;
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys, 
    const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
	vec3 val = 0.0;
	double u = 0.0;

	curve.clear();

	int numSegments = keys.size() - 1;
	for (int segment = 0; segment < numSegments; segment++)
    {
        for (double t = keys[segment].first; t < keys[segment+1].first - FLT_EPSILON; t += mDt)
        {
			// TODO: Compute u, fraction of duration between segment and segment next, for example,
			// u = 0.0 when t = keys[segment-1].first  
			// u = 1.0 when t = keys[segment].first

			//if (segment > 0) {
			//
			//	u = (t - keys[segment - 1].first) / (keys[segment].first - keys[segment - 1].first);
			//}
			u = (t - keys[segment].first) / (keys[segment + 1].first - keys[segment].first);

            val = interpolateSegment(keys, ctrlPoints, segment, u);
            curve.push_back(val);
        }
    }
	// add last point
	if (keys.size() > 1)
	{
		u = 1.0;
		val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
		curve.push_back(val);
	}
}


// Interpolate p0 and p1 so that t = 0 returns p0 and t = 1 returns p1
vec3 ALinearInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO: Linear interpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1

	curveValue = key0 * (1 - u) + key1 * u;


	return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	curveValue = pow((1 - u), 3) * b0 + 3 * u * pow((1 - u), 2) * b1 + 3 * pow(u, 2) * b2 * (1 - u) + pow(u, 3) * b3;

	return curveValue;
}

vec3 ACasteljauInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	vec3 b_1_0 = LerpVec3(b0, b1, u);
	vec3 b_1_1 = LerpVec3(b1, b2, u);
	vec3 b_1_2 = LerpVec3(b2, b3, u);
	vec3 b_2_0 = LerpVec3(b_1_0, b_1_1, u);
	vec3 b_2_1 = LerpVec3(b_1_1, b_1_2, u);
	vec3 b_3_0 = LerpVec3(b_2_0, b_2_1, u);

	curveValue = b_3_0;
	
	return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
	// Hint: Using Eigen::MatrixXd data representations for a matrix operations

	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];

	Eigen::MatrixXd M(4,4);
	M(0, 0) = 1; M(0, 1) = -3; M(0, 2) = 3; M(0, 3) = -1; 
	M(1, 0) = 0; M(1, 1) = 3; M(1, 2) = -6; M(1, 3) = 3;
	M(2, 0) = 0; M(2, 1) = 0; M(2, 2) = 3; M(2, 3) = -3;
	M(3, 0) = 0; M(3, 1) = 0; M(3, 2) = 0; M(3, 3) = 1;

	Eigen::MatrixXd G(3, 4);
	G(0, 0) = b0[0]; G(0, 1) = b1[0]; G(0, 2) = b2[0]; G(0, 3) = b3[0];
	G(1, 0) = b0[1]; G(1, 1) = b1[1]; G(1, 2) = b2[1]; G(1, 3) = b3[1];
	G(2, 0) = b0[2]; G(2, 1) = b1[2]; G(2, 2) = b2[2]; G(2, 3) = b3[2];

	Eigen::MatrixXd V(3, 4);
	V = G * M;

	Eigen::MatrixXd U(4, 1);
	U(0, 0) = 1;
	U(1, 0) = u;
	U(2, 0) = pow(u, 2);
	U(3, 0) = pow(u, 3);

	Eigen::MatrixXd fu(3, 1);
	fu = V * U;

	curveValue[0] = fu(0, 0);
	curveValue[1] = fu(1, 0);
	curveValue[2] = fu(2, 0);

	return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 p0;
	vec3 p1;
	vec3 q0; // slope at p0
	vec3 q1; // slope at p1
	vec3 curveValue(0, 0, 0);

	// TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  
	p0 = keys[segment].second;
	p1 = keys[segment + 1].second;
	q0 = ctrlPoints[segment];
	q1 = ctrlPoints[segment + 1];

	curveValue = (2 * pow(u, 3) - 3 * pow(u, 2) + 1) * p0
		+ (-2 * pow(u, 3) + 3 * pow(u, 2)) * p1
		+ (pow(u, 3) - 2 * pow(u, 2) + u) * q0
		+ (pow(u, 3) - pow(u, 2)) * q1;


	return curveValue;
}

vec3 ABSplineInterpolatorVec3::interpolateSegment(
    const std::vector<ASplineVec3::Key>& keys,
    const std::vector<vec3>& ctrlPoints, 
    int segment, double u)
{
	vec3 curveValue(0, 0, 0);

	// Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = curve interval on knot vector in which to interpolate
	//     t = time value	

	// Step 1: determine the index j
	// Step 2: compute the n nonzero Bspline Basis functions N given j
	// Step 3: get the corresponding control points from the ctrlPoints vector
	// Step 4: compute the Bspline curveValue at time t

	return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys, 
    std::vector<vec3>& ctrlPoints, 
    vec3& startPoint, vec3& endPoint)
{
    ctrlPoints.clear();
    if (keys.size() <= 1) return;

    for (int i = 1; i < keys.size(); i++)
    {
        vec3 b0, b1, b2, b3;
		// TODO: compute b0, b1, b2, b3

		// when i = 1, start at segment 0
		if (i < (keys.size() - 1) && i > 1) {
		vec3 p0 = keys[i - 2].second;
		vec3 p1 = keys[i - 1].second;
		vec3 p2 = keys[i].second;
		vec3 p3 = keys[i + 1].second;
		double t0 = keys[i - 2].first;
		double t1 = keys[i - 1].first;
		double t2 = keys[i].first;
		double t3 = keys[i + 1].first;

		// slope at p1 by averaging slope on the left and right
		vec3 s_left_p1 = (p1 - p0) / (t1 - t0);
		vec3 s_right_p1 = (p2 - p1) / (t2 - t1);
		vec3 s1 = (s_left_p1 + s_right_p1) / 2;

		// slope at p2 by averaging slope on the left and right
		vec3 s_left_p2 = (p2 - p1) / (t2 - t1);
		vec3 s_right_p2 = (p3 - p2) / (t3 - t2);
		vec3 s2 = (s_left_p2 + s_right_p2) / 2;

		// find b
		b0 = p1;
		b3 = p2;

		b1 = b0 + s1 / 3;
		b2 = b3 - s2 / 3;
		
		} // left end point with more than one segment
		else if (i == 1 && i != (keys.size() - 1)) {
			vec3 p1 = keys[i - 1].second;
			vec3 p2 = keys[i].second;
			vec3 p3 = keys[i + 1].second;

			double t1 = keys[i - 1].first;
			double t2 = keys[i].first;
			double t3 = keys[i + 1].first;

			// slope at p1 is only the slope from right
			vec3 s_right_p1 = (p2 - p1) / (t2 - t1);
			vec3 s1 = (s_right_p1) / 1;

			// slope at p2 by averaging slope on the left and right
			vec3 s_left_p2 = (p2 - p1) / (t2 - t1);
			vec3 s_right_p2 = (p3 - p2) / (t3 - t2);
			vec3 s2 = (s_left_p2 + s_right_p2) / 2;

			// find b
			b0 = p1;
			b3 = p2;

			b1 = b0 + s1 / 3;
			b2 = b3 - s2 / 3;

		} // right end point with more than one segment
		else if (i == (keys.size() - 1) && i != 1) {
			vec3 p0 = keys[i - 2].second;
			vec3 p1 = keys[i - 1].second;
			vec3 p2 = keys[i].second;
			double t0 = keys[i - 2].first;
			double t1 = keys[i - 1].first;
			double t2 = keys[i].first;

			// slope at p1 by averaging slope on the left and right
			vec3 s_left_p1 = (p1 - p0) / (t1 - t0);
			vec3 s_right_p1 = (p2 - p1) / (t2 - t1);
			vec3 s1 = (s_left_p1 + s_right_p1) / 2;

			// slope at p2 only the slope from left
			vec3 s_left_p2 = (p2 - p1) / (t2 - t1);
			vec3 s2 = (s_left_p2 ) / 1;

			// find b
			b0 = p1;
			b3 = p2;

			b1 = b0 + s1 / 3;
			b2 = b3 - s2 / 3;
		}// only one segment
		else if (i == (keys.size() - 1) && i == 1) {
			vec3 p1 = keys[i - 1].second;
			vec3 p2 = keys[i].second;
			double t1 = keys[i - 1].first;
			double t2 = keys[i].first;

			// slope at p1 by averaging slope on the left and right
			vec3 s_right_p1 = (p2 - p1) / (t2 - t1);
			vec3 s1 = (s_right_p1) / 1;

			// slope at p2 only the slope from left
			vec3 s_left_p2 = (p2 - p1) / (t2 - t1);
			vec3 s2 = (s_left_p2) / 1;

			// find b
			b0 = p1;
			b3 = p2;

			b1 = b0 + s1 / 3;
			b2 = b3 - s2 / 3;
		}
		


        ctrlPoints.push_back(b0);
        ctrlPoints.push_back(b1);
        ctrlPoints.push_back(b2);
        ctrlPoints.push_back(b3);
    }
}

void AHermiteInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys,
    std::vector<vec3>& ctrlPoints,
    vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	ctrlPoints.resize(keys.size(), vec3(0, 0, 0));
	if (keys.size() <= 1) return;

	// TODO: 
	// For each key point pi, compute the corresonding value of the slope pi_prime.
	// Hints: Using Eigen::MatrixXd for a matrix data structures, 
	// this can be accomplished by solving the system of equations AC=D for C.
	// Don't forget to save the values computed for C in ctrlPoints
	// For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
	// For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

	// Step 1: Initialize A
	// Step 2: Initialize D
	// Step 3: Solve AC=D for C
	// Step 4: Save control points in ctrlPoints

	// Control Points: [p0_prime, p1_prime, p2_prime, ..., pm_prime]
	// The size of control points should be the same as the size of the keys
	// Use operator[] to set elements in ctrlPoints by indices

	// Hint: Do not use push_back() to insert control points here because the vector has been resized
	

	// step 0 choose natural or clamp end points;
	//bool ifNatural = true;
	//Step 1: Initialize A
	int N = keys.size() - 1;
	Eigen::MatrixXd A(N + 1, N + 1);

	//clamp end points
	if (!ifNatural) {
		for (int row = 0; row < N + 1; row++) {
			for (int col = 0; col < N + 1; col++) {

				// first row of A
				if (row == 0) {
					if (col == 0) {
						A(row, col) = 1;
					}
					else {
						A(row, col) = 0;
					}
				}// last row of A
				else if (row == N) {
					if (col == N) {
						A(row, col) = 1;
					}
					else {
						A(row, col) = 0;
					}
				}// other rows
				else {
					// insert first 1 in each row
					if (col == (row - 1)) {
						A(row, col) = 1;
					}
					else if (col == (row)) {
						A(row, col) = 4;
					}
					else if (col == (row + 1)) {
						A(row, col) = 1;
					}
					else {
						A(row, col) = 0;
					}
				}
			}
		}
	}
	// natural end points
	else {
		for (int row = 0; row < N + 1; row++) {
			for (int col = 0; col < N + 1; col++) {
				// first row of A
				if (row == 0) {
					if (col == 0) {
						A(row, col) = 4;
					}
					else if (col == 1) {
						A(row, col) = 2;
					}
					else {
						A(row, col) = 0;
					}
				}// last row of A
				else if (row == N) {
					if (col == N - 1) {
						A(row, col) = 2;
					}
					else if (col == N) {
						A(row, col) = 4;
					}
					else {
						A(row, col) = 0;
					}
				}// other rows
				else {
					// insert first 1 in each row
					if (col == (row - 1)) {
						A(row, col) = 1;
					}
					else if (col == (row)) {
						A(row, col) = 4;
					}
					else if (col == (row + 1)) {
						A(row, col) = 1;
					}
					else {
						A(row, col) = 0;
					}
				}
			}
		}
	}

	// Step 2: Initialize D
	//vec3 s0 = (keys[1].second - keys[0].second) / (keys[1].first - keys[0].first);
	//vec3 s1 = (keys[N].second - keys[N - 1].second) / (keys[N].first - keys[N - 1].first);
	vec3 s0 = (keys[1].second - keys[0].second) / 1;
	vec3 s1 = (keys[N].second - keys[N - 1].second) / 1;

	Eigen::MatrixXd D(N + 1, 3);

	// clamp end points
	if (!ifNatural) {
		for (int row = 0; row < N + 1; row++) {
			// first row
			if (row == 0) {
				D(row, 0) = s0[0];
				D(row, 1) = s0[1];
				D(row, 2) = s0[2];
			}// last row of A
			else if (row == N) {
				D(row, 0) = s1[0];
				D(row, 1) = s1[1];
				D(row, 2) = s1[2];
			}// other rows
			else {
				D(row, 0) = (3 * (keys[row + 1].second - keys[row - 1].second))[0];
				D(row, 1) = (3 * (keys[row + 1].second - keys[row - 1].second))[1];
				D(row, 2) = (3 * (keys[row + 1].second - keys[row - 1].second))[2];
			}
		}
	}
	// natrual end points
	else {
		for (int row = 0; row < N + 1; row++) {
			// first row
			if (row == 0) {
				D(row, 0) = (6 * (keys[row + 1].second - keys[row].second))[0];
				D(row, 1) = (6 * (keys[row + 1].second - keys[row].second))[1];
				D(row, 2) = (6 * (keys[row + 1].second - keys[row].second))[2];
			}// last row of A
			else if (row == N) {
				D(row, 0) = (6 * (keys[row].second - keys[row - 1].second))[0];
				D(row, 1) = (6 * (keys[row].second - keys[row - 1].second))[1];
				D(row, 2) = (6 * (keys[row].second - keys[row - 1].second))[2];
			}// other rows
			else {
				D(row, 0) = (3 * (keys[row + 1].second - keys[row - 1].second))[0];
				D(row, 1) = (3 * (keys[row + 1].second - keys[row - 1].second))[1];
				D(row, 2) = (3 * (keys[row + 1].second - keys[row - 1].second))[2];
			}
		}
	}


	// Step 3: Solve AC=D for C
	Eigen::MatrixXd C(N + 1, 3);

	C = (A.inverse()) * D;

	// Step 4: Save control points in ctrlPoints
	for (int i = 0; i < keys.size(); i++) {
		ctrlPoints[i][0] = C(i, 0);
		ctrlPoints[i][1] = C(i, 1);
		ctrlPoints[i][2] = C(i, 2);
	}


}

void ABSplineInterpolatorVec3::computeControlPoints(
    const std::vector<ASplineVec3::Key>& keys,
    std::vector<vec3>& ctrlPoints, 
    vec3& startPt, vec3& endPt)
{
    ctrlPoints.clear();
	ctrlPoints.resize(keys.size() + 2, vec3(0, 0, 0));
    if (keys.size() <= 1) return;

	// TODO:
	// Hints: 
	// 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

	// 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = interval on knot vector in which to interpolate
	//     t = time value
	//     l = derivative (l = 1 => 1st derivative)

	// Step 1: Calculate knot vector using a uniform BSpline
	//         (assune knots are evenly spaced 1 apart and the start knot is at time = 0.0)

	// Step 2: Calculate A matrix  for a natural BSpline
	//         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)

	// Step 3: Calculate  D matrix composed of our target points to interpolate

	// Step 4: Solve AC=D for C 

	// Step 5: save control points in ctrlPoints

	// Hint: Do not use push_back() to insert control points here because the vector has been resized



	//

}


vec3 AEulerLinearInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys, 
	const std::vector<vec3>& ctrlPoints, 
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO:
	// Linear interpolate between key0 and key1
	// You should convert the angles to find the shortest path for interpolation

	return curveValue;
}

vec3 AEulerCubicInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys, 
	const std::vector<vec3>& ctrlPoints, int segment, double t)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
	// You should convert the angles to find the shortest path for interpolation

	return curveValue;
}

void AEulerCubicInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys, 
	std::vector<vec3>& ctrlPoints, vec3 & startPoint, vec3 & endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	// Hint: One naive way is to first convert the keys such that the differences of the x, y, z Euler angles 
	//		 between every two adjacent keys are less than 180 degrees respectively 

	for (int i = 1; i < keys.size(); i++)
	{
		vec3 b0, b1, b2, b3;

		// TODO: compute b0, b1, b2, b3

		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
	}
}


