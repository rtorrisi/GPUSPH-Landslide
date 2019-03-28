#ifndef _CHUTELANDSLIDE_H
#define	_CHUTELANDSLIDE_H

#include "XProblem.h"

#define __toRadians(radius) radius*M_PI/180.f

class ChuteLandslide: public XProblem {
    private:
		/*
			All measures are in meters. Angles are specified in degree, then converted in radiants;

		*/
		double		m_sphereRadius = 0.195f;
		double		m_cupLatitudeCutAngle = __toRadians(30);
		double		m_cupInitialHeight = m_sphereRadius - m_sphereRadius * sin(m_cupLatitudeCutAngle);
		double		m_chuteInclinationAngle = __toRadians(45);
		double		m_chuteWidth = 0.8f; //1.6f;
		double		m_chuteThickness = 0.f; //init after set_deltap()
		double		m_chuteHorizontalLength = 1.f; //2.25f;
		double		m_chuteHorizontalOffset = -0.04; //-0.095f;
		double		m_chuteObliqueLength = 1.56f;
		double		m_chuteObliqueDelta = 0.f; //init after set_deltap()

    public:
        ChuteLandslide(GlobalData *);
		void setFramework();
		void setSPHParameters();
		void setPhysicalParameters();
		void buildGeometry();
};

#endif /*_CHUTELANDSLIDE_H */