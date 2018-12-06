#ifndef _LANDSLIDE_H
#define	_LANDSLIDE_H

#include "XProblem.h"

class Landslide: public XProblem {
    private:
		double		m_chuteInclinationAngle = (5*M_PI)/4.f;
		double		m_fluidRadius = 0.195;
		double		m_chuteThickness = 0.2; //0.0035
		double		m_chuteWidth = 1.6;
		double		m_chuteHorizontalLength = 2.25;
		double		m_chuteObliqueLength = 1.56;
		double		m_deltaChuteObliqueLength = sqrt(2)*m_chuteThickness;
		double		m_offset = -0.095;

    public:
        Landslide(GlobalData *);
		void setFramework();
		void setSPHParameters();
		void setPhysicalParameters();
		void buildGeometry();
};

#endif /*_LANDSLIDE_H */