 #include "ASplineQuat.h"
#include <algorithm>
#pragma warning(disable:4018)

ASplineQuat::ASplineQuat() : mDt(1.0 / 120.0), mLooping(true), mType(LINEAR)
{
}

ASplineQuat::~ASplineQuat()
{
}

void ASplineQuat::setInterpolationType(ASplineQuat::InterpolationType type)
{
    mType = type;
    cacheCurve();
}

ASplineQuat::InterpolationType ASplineQuat::getInterpolationType() const
{
    return mType;
}

void ASplineQuat::setLooping(bool loop)
{
    mLooping = loop;
}

bool ASplineQuat::getLooping() const
{
    return mLooping;
}

void ASplineQuat::setFramerate(double fps)
{
    mDt = 1.0 / fps;
}

double ASplineQuat::getFramerate() const
{
    return 1.0 / mDt;
}

int ASplineQuat::getCurveSegment(double time)
{
	int segment = 0;
	bool foundSegment = false;

	double t = time;
	if (t < 0.0)
		t = 0.0;

	int numKeys = mKeys.size();
	while (!foundSegment) {
		if (segment == numKeys - 1) {
			segment = numKeys - 2;
			foundSegment = true;
		}
		else {
			double keyTime0 = mKeys[segment].first;
			double keyTime1 = mKeys[segment + 1].first;
			if ((t >= keyTime0) && (t < keyTime1))
				 foundSegment = true;
			else segment++;
		}
	}
	return segment;

}


quat ASplineQuat::getCachedValue(double t) const
{

	if (mCachedCurve.empty() || mKeys.empty()) return quat();

	if (t < mKeys[0].first)
		return mCachedCurve[0];
	else
		t -= mKeys[0].first;

	int numFrames = (int)(t / mDt);
	int i = mLooping ? numFrames % mCachedCurve.size() : std::min<int>(numFrames, mCachedCurve.size() - 1);
	int inext = mLooping ? (i + 1) % mCachedCurve.size() : std::min<int>(i + 1, mCachedCurve.size() - 1);
	quat key1 = mCachedCurve[i];
	quat key2 = mCachedCurve[inext];
	double u = (t - numFrames * mDt) / mDt;
	return quat::Slerp(key1, key2, u);

}

void ASplineQuat::cacheCurve()
{
	int numKeys = mKeys.size();

	if (numKeys == 1)
	{
		mCachedCurve.clear();
		mCachedCurve.push_back(mKeys[0].second);
	}

	if (mType == LINEAR && numKeys >= 2)
		createSplineCurveLinear();

	if (mType == CUBIC && numKeys >= 2)
	{
		quat startQuat = mKeys[0].second;
		quat endQuat = mKeys[numKeys-1].second;

		computeControlPoints(startQuat, endQuat);
		createSplineCurveCubic();
	}
}
void ASplineQuat::computeControlPoints(quat& startQuat, quat& endQuat)
{
	// startQuat is a phantom point at the left-most side of the spline
	// endQuat is a phantom point at the left-most side of the spline

	mCtrlPoints.clear();
	int numKeys = mKeys.size();
	if (numKeys <= 1) return;

	quat b0, b1, b2, b3;
	quat q_1, q0, q1, q2;

	for (int segment = 0; segment < numKeys - 1; segment++)
	{
		// TODO: student implementation goes here
		//  Given the quaternion keys q_1, q0, q1 and q2 associated with a curve segment, compute b0, b1, b2, b3 
		//  for each cubic quaternion curve, then store the results in mCntrlPoints in same the same way 
		//  as was used with the SplineVec implementation
		//  Hint: use the SDouble, SBisect and Slerp to compute b1 and b2

		// middle segment with more than one segment
		std::cout << "===== segment: " << segment <<" =======" << std::endl;

		if (segment > 0 && segment != numKeys - 2) {
			q_1 = mKeys[segment - 1].second;
			q0 = mKeys[segment].second;
			q1 = mKeys[segment + 1].second;
			q2 = mKeys[segment + 2].second;

			quat q_prim0 = quat::SDouble(q2, q1);
			quat q_prim1 = quat::SDouble(q_1, q0);

			quat q_star0 = quat::SBisect(q0, q_prim0);
			quat q_star1 = quat::SBisect(q_prim1, q1);

			b0 = q0;
			b1 = quat::Slerp(q0, q_star1, 1.0 / 3.0);
			b2 = quat::Slerp(q1, q_star0, 1.0 / 3.0);
			b3 = q1;
			std::cout << "middle" << std::endl;

			//std::cout << "q_star0: " << q_star0[3] << std::endl;
			//std::cout << "q_star0: " << q_star0[0] << std::endl;
			//std::cout << "q_star0: " << q_star0[1] << std::endl;
			//std::cout << "q_star0: " << q_star0[2] << std::endl;
			//std::cout << "q_star1: " << q_star1[3] << std::endl;
			//std::cout << "q_star1: " << q_star1[0] << std::endl;
			//std::cout << "q_star1: " << q_star1[1] << std::endl;
			//std::cout << "q_star1: " << q_star1[2] << std::endl;


			
		}// left end point with more than one segment
		else if(segment == 0 && segment != numKeys - 2) {
			
			q0 = mKeys[segment].second;
			q1 = mKeys[segment + 1].second;
			q2 = mKeys[segment + 2].second;
			q_1 = quat::SDouble(q1, q0);

			quat q_prim0 = quat::SDouble(q2, q1);
			quat q_prim1 = quat::SDouble(q_1, q0);

			quat q_star0 = quat::SBisect(q0, q_prim0);
			quat q_star1 = quat::SBisect(q_prim1, q1);

			b0 = q0;
			b1 = quat::Slerp(q0, q_star1, 1.0 / 3.0);
			b2 = quat::Slerp(q1, q_star0, 1.0 / 3.0);
			b3 = q1;
			std::cout << "left" << std::endl;
		
		}// right end point with more than one segment
		else if (segment == numKeys - 2 && segment != 0) {
			q_1 = mKeys[segment - 1].second;
			q0 = mKeys[segment].second;
			q1 = mKeys[segment + 1].second;
			q2 = quat::SDouble(q0, q1);

			quat q_prim0 = quat::SDouble(q2, q1);
			quat q_prim1 = quat::SDouble(q_1, q0);

			quat q_star0 = quat::SBisect(q0, q_prim0);
			quat q_star1 = quat::SBisect(q_prim1, q1);

			b0 = q0;
			b1 = quat::Slerp(q0, q_star1, 1.0 / 3.0);
			b2 = quat::Slerp(q1, q_star0, 1.0 / 3.0);
			b3 = q1;
			std::cout << "right" << std::endl;

		}// only one segment
		else if (segment == numKeys - 2 && segment == 0) {
			q0 = mKeys[segment].second;
			q1 = mKeys[segment + 1].second;
			q_1 = quat::SDouble(q1, q0);
			q2 = quat::SDouble(q0, q1);

			quat q_prim0 = quat::SDouble(q2, q1);
			quat q_prim1 = quat::SDouble(q_1, q0);

			quat q_star0 = quat::SBisect(q0, q_prim0);
			quat q_star1 = quat::SBisect(q_prim1, q1);

			b0 = q0;
			b1 = quat::Slerp(q0, q_star1, 1.0 / 3.0);
			b2 = quat::Slerp(q1, q_star0, 1.0 / 3.0);
			b3 = q1;
			//std::cout << "calculated" << std::endl;

		}
		std::cout << "b0VW: " << b0[3]<< std::endl;
		std::cout << "b0VX: " << b0[0] << std::endl;
		std::cout << "b0VY: " << b0[1] << std::endl;
		std::cout << "b0VZ: " << b0[2] << std::endl;
		std::cout << "b1VW: " << b1[3] << std::endl;
		std::cout << "b1VX: " << b1[0] << std::endl;
		std::cout << "b1VY: " << b1[1] << std::endl;
		std::cout << "b1VZ: " << b1[2] << std::endl;
		std::cout << "b2VW: " << b2[3] << std::endl;
		std::cout << "b2VX: " << b2[0] << std::endl;
		std::cout << "b2VY: " << b2[1] << std::endl;
		std::cout << "b2VZ: " << b2[2] << std::endl;
		std::cout << "b3VW: " << b3[3] << std::endl;
		std::cout << "b3VX: " << b3[0] << std::endl;
		std::cout << "b3VY: " << b3[1] << std::endl;
		std::cout << "b3VZ: " << b3[2] << std::endl;
		



		mCtrlPoints.push_back(b0);
		mCtrlPoints.push_back(b1);
		mCtrlPoints.push_back(b2);
		mCtrlPoints.push_back(b3);
	}
}

quat ASplineQuat::getLinearValue(double t)
{

	quat q;
	int segment = getCurveSegment(t);

	// TODO: student implementation goes here
	// compute the value of a linear quaternion spline at the value of t using slerp
	quat q0 = mKeys[segment].second;
	quat q1 = mKeys[segment + 1].second;


	double t0 = mKeys[segment].first;
	double t1 = mKeys[segment + 1].first;

	double u = (t - t0) / (t1 - t0);

	q = quat::Slerp(q0, q1, u);

	//std::cout << "segment: "<<segment << std::endl;
	//std::cout << "mKeys size: " << mKeys.size() << std::endl;
	//
	//std::cout << "VW: " << q[3] << std::endl;
	//std::cout << "VX: " << q[0] << std::endl;
	//std::cout << "VY: " << q[1] << std::endl;
	//std::cout << "VZ: " << q[2] << std::endl;
	
	return q;	
}

void ASplineQuat::createSplineCurveLinear()
{
	std::cout << "createSplineCurveLinear" << std::endl;
	quat q;
	mCachedCurve.clear();
	int numKeys = mKeys.size(); 
	double startTime = mKeys[0].first;
	double endTime = mKeys[numKeys-1].first;

	for (double t = startTime; t <= endTime; t += mDt)
	{
		q = getLinearValue(t);
		mCachedCurve.push_back(q);
	}

	for (int i = 0; i < mKeys.size(); i++) {
	
		std::cout << "time: " << mKeys[i].first << " qVW: " << mKeys[i].second[3] << " qVX: " << mKeys[i].second[0] << " qVY: " << mKeys[i].second[1] << " qVW: " << mKeys[i].second[2] << std::endl;
	}
}


quat ASplineQuat::getCubicValue(double t)
{
	quat q, b0, b1, b2, b3;
	int segment = getCurveSegment(t);

	// TODO: student implementation goes here
	// compute the value of a cubic quaternion spline at the value of t using Scubic
	b0 = mCtrlPoints[4 * segment];
	b1 = mCtrlPoints[4 * segment + 1];
	b2 = mCtrlPoints[4 * segment + 2];
	b3 = mCtrlPoints[4 * segment + 3];
	double t0 = mKeys[segment].first;
	double t1 = mKeys[segment + 1].first;

	double u = (t - t0) / (t1 - t0);

	q = quat::Scubic(b0, b1, b2, b3, u);

	//std::cout << "====================="<< std::endl;
	//std::cout << "q: " << q[3] << std::endl;
	//std::cout << "q: " << q[0] << std::endl;
	//std::cout << "q: " << q[1] << std::endl;
	//std::cout << "q: " << q[2] << std::endl;

	return q;
}

void ASplineQuat::createSplineCurveCubic()
{
	quat q;
	mCachedCurve.clear();
	int numKeys = mKeys.size();
	double startTime = mKeys[0].first;
	double endTime = mKeys[numKeys - 1].first;

	for (double t = startTime; t <= endTime; t += mDt)
	{
		q = getCubicValue(t);
		mCachedCurve.push_back(q);
	}
}


void ASplineQuat::editKey(int keyID, const quat& value)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys[keyID].second = value;
	cacheCurve();
}

void ASplineQuat::appendKey(const quat& value, bool updateCurve)
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
	for (int i = 0; i < mKeys.size(); i++) {

		//std::cout << "time: " << mKeys[i].first << " qVW: " << mKeys[i].second[3]
		//	<< " qVX: " << mKeys[i].second[0] << " qVY: " << mKeys[i].second[1] << " qVZ: " << mKeys[i].second[2] << std::endl;
	}
}

int ASplineQuat::insertKey(double time, const quat& value, bool updateCurve)
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
			if (updateCurve) cacheCurve();
			return i;
		}
	}
	// Append at the end of the curve
	appendKey(time, value, updateCurve);
	return mKeys.size() - 1;
}

void ASplineQuat::appendKey(double t, const quat& value, bool updateCurve)
{
    mKeys.push_back(Key(t, value));
    if (updateCurve) cacheCurve();
}

void ASplineQuat::deleteKey(int keyID)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    mKeys.erase(mKeys.begin() + keyID);
	cacheCurve();
}

quat ASplineQuat::getKey(int keyID)
{
    assert(keyID >= 0 && keyID < mKeys.size());
    return mKeys[keyID].second;
}

int ASplineQuat::getNumKeys() const
{
    return mKeys.size();
}

void ASplineQuat::clear()
{
    mKeys.clear();
}

double ASplineQuat::getDuration() const
{
    return mCachedCurve.size() * mDt;
}

double ASplineQuat::getNormalizedTime(double t) const
{
    double duration = getDuration();
    int rawi = (int)(t / duration);
    return t - rawi*duration;
}
