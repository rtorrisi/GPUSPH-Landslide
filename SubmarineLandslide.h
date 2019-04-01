#ifndef _SUBMARINELANDSLIDE_H
#define	_SUBMARINELANDSLIDE_H

#include "XProblem.h"

#define __toRadians(radius) radius*M_PI/180.f

class SubmarineLandslide: public XProblem {
    private:
		/*
			All measures are in meters. Angles are specified in degree, then converted in radiants;
		*/

		double		m_waterBoxWidth			= 4.f;
		double		m_waterBoxHeight		= 1.6f;
		
		double		m_bulkWidth				= 0.65f;
		double		m_bulkHeight			= 0.65f;

		double		m_chuteInclinationAngle = __toRadians(45);
		double		m_chuteUpperWidth		= 1.f;
		double		m_chuteHeight			= 1.5f;
		double		m_chuteWidth 			= 2.5f;
		double		m_chuteThickness		= 0.15f;

		double		m_depth					= 1.f;
		double		m_gapFixOffset			= 0.025f;
		double		m_gapFixOffset2			= 0.005f;

    public:
        SubmarineLandslide(GlobalData *);
		void setFramework();
		void setSPHParameters();
		void setPhysicalParameters();
		void buildGeometry();
};

#endif /*_SUBMARINELANDSLIDE_H */