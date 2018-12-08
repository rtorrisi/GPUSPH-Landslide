#include "Landslide.h"
#include "GlobalData.h"
#include "cudasimframework.cu"

Landslide::Landslide(GlobalData *_gdata) : XProblem(_gdata)
{
	m_name = "Landslide";

	setFramework();

	const int mlsIters = get_option("mls", 0); // --mls N to enable MLS filter every N iterations
	const int ppH = get_option("ppH", 16); // --ppH N to change deltap to H/N

	if (mlsIters > 0) addFilter(MLS_FILTER, mlsIters);

	set_deltap(m_cupInitialHeight/ppH);
	
	m_chuteThickness = 3*m_deltap;
	m_chuteObliqueDelta = sqrt(2) * m_chuteThickness;
	m_chuteObliqueLength += m_chuteObliqueDelta;

	// SPH parameters
	setSPHParameters();

	// Physical parameters
	setPhysicalParameters();

	// Drawing and saving times
	add_writer(VTKWRITER, 0.01);

	// Building the geometry
	buildGeometry();
}

void Landslide::setFramework()
{
	// density diffusion terms, see DensityDiffusionType
	const int	rhodiff = get_option("density-diffusion", 1);

	SETUP_FRAMEWORK(
		viscosity<DYNAMICVISC>,
		boundary<DYN_BOUNDARY>,
		add_flags<ENABLE_PLANES>
	);
}

void Landslide::setSPHParameters()
{
	simparams()->buildneibsfreq = 10;
}

void Landslide::setPhysicalParameters()
{
	physparams()->gravity = make_float3(0.0, 0.0, -9.81f);
	const float g = length(physparams()->gravity);
	const float maxvel = sqrt(2*g*1);
	// purely for cosmetic reason, let's round the soundspeed to the next integer
	const float c0 = ceil(10*maxvel);
	add_fluid(1000.0);
	set_equation_of_state(0, 7.0f, c0);
	physparams()->dcoeff = 5.0f*g*1;
	physparams()->r0 = m_deltap;
	//physparams()->visccoeff = 0.05f;
	set_kinematic_visc(0, 3.0e-2f);
	//set_kinematic_visc(0, 1.0e-6f);
	physparams()->artvisccoeff = 0.3f;
	physparams()->epsartvisc = 0.01*simparams()->slength*simparams()->slength;
	physparams()->epsxsph = 0.5f;
}

void Landslide::buildGeometry()
{
	setPositioning(PP_CENTER);

	double3 maxChutePoint = make_double3( 0, //x
		- m_chuteObliqueLength*cos(m_chuteInclinationAngle), //y
		m_chuteObliqueLength*sin(m_chuteInclinationAngle) //z
	);

	double3 cupPositionDelta = make_double3( 0, //x
		(m_sphereRadius*cos(m_cupLatitudeCutAngle))*cos(m_chuteInclinationAngle), //y
		- (m_sphereRadius*cos(m_cupLatitudeCutAngle))*sin(m_chuteInclinationAngle) //z
	);

	double3 cupLatitudeCutDelta = make_double3( 0, //x
		- (m_sphereRadius-m_cupInitialHeight)*cos(m_chuteInclinationAngle), //y
		- (m_sphereRadius-m_cupInitialHeight)*sin(m_chuteInclinationAngle) //z
	);

	Point cupOriginPoint = Point(
		m_chuteWidth / 2.f,
		maxChutePoint.y + cupPositionDelta.y + cupLatitudeCutDelta.y,
		maxChutePoint.z + cupPositionDelta.z + cupLatitudeCutDelta.z
	);

	GeometryID fluid = addSphere(GT_FLUID, FT_SOLID, cupOriginPoint, m_sphereRadius);
	rotate(fluid, m_chuteInclinationAngle, 0.f, 0.f);

	GeometryID obliquePlane = addPlane(0.0f, 1.0f, 1.0f, 0.0f);
	setIntersectionType(obliquePlane, IT_INTERSECT);
	//deleteGeometry(obliquePlane);

	setPositioning(PP_CORNER);

	GeometryID horizontalBox = addBox(GT_FIXED_BOUNDARY, FT_SOLID,
		Point(0, m_chuteHorizontalOffset, 0), m_chuteWidth, m_chuteHorizontalLength, m_chuteThickness);

	GeometryID obliqueBox = addBox(GT_FIXED_BOUNDARY, FT_SOLID,
		Point(0,0,0), m_chuteWidth, m_chuteObliqueLength, m_chuteThickness);
	rotate(obliqueBox, M_PI+m_chuteInclinationAngle, 0.f, 0.f);

	GeometryID horizontalPlane = addPlane(0.f, 0.f, 1.f, 0.01f);
	setIntersectionType(horizontalPlane, IT_INTERSECT);
	//deleteGeometry(horizontalPlane);
}