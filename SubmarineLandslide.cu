#include <string>

#include "SubmarineLandslide.h"
#include "GlobalData.h"
#include "cudasimframework.cu"

SubmarineLandslide::SubmarineLandslide(GlobalData *_gdata) : XProblem(_gdata)
{
	m_name = "SubmarineLandslide";

	setFramework();

	const int mlsIters = get_option("mls", 0); // --mls N to enable MLS filter every N iterations
	const int ppH = get_option("ppH", 32); // --ppH N to change deltap to H/N

	if (mlsIters > 0) addFilter(MLS_FILTER, mlsIters);

	set_deltap(m_bulkHeight/ppH);

	// SPH parameters
	setSPHParameters();

	// Physical parameters
	setPhysicalParameters();

	// Drawing and saving times
	add_writer(VTKWRITER, 0.01);

	// Building the geometry
	buildGeometry();
}

void SubmarineLandslide::setFramework()
{
	// density diffusion terms, see DensityDiffusionType
	const int	rhodiff = get_option("density-diffusion", 1);
	const bool	newtonian = get_option("newtonian", false);

	SETUP_FRAMEWORK(
		rheology<PAPANASTASIOU>,
		turbulence_model<LAMINAR_FLOW>,
		visc_model<MORRIS>,
		visc_average<HARMONIC>,
		computational_visc<DYNAMIC>,
		boundary<DYN_BOUNDARY>,
		periodicity<PERIODIC_Y>
	).select_options(
		newtonian, rheology<NEWTONIAN>()
	);
}

void SubmarineLandslide::setSPHParameters()
{
	simparams()->buildneibsfreq = 10;
	simparams()->tend = 5.0f;
}

void SubmarineLandslide::setPhysicalParameters()
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
	if (YIELDING_RHEOLOGY(simparams()->rheologytype))
		set_yield_strength(0, 1.0f);
}

void SubmarineLandslide::buildGeometry()
{
	setPositioning(PP_CORNER);

	GeometryID obliqueChute = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(0.f, 0.f, 0.f),
		m_chuteWidth, m_depth, m_chuteHeight);

	// This plane is used to erase unnecessary fixed boundary particles
	GeometryID optimizationPlane = addPlane(1.f, 0.f, 1.f, -m_chuteWidth + m_chuteThickness, FT_UNFILL);
		setIntersectionType(optimizationPlane, IT_INTERSECT);

	GeometryID horizontalChute = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(0.f, 0.f, 0.f),
		m_chuteUpperWidth - m_gapFixOffset2, m_depth, m_chuteHeight);

	// This box is used to erase unnecessary fixed boundary particles
	GeometryID optimizationBox = addBox(GT_FIXED_BOUNDARY, FT_NOFILL, Point(0.f, 0.f, 0.f),
		m_chuteUpperWidth - m_gapFixOffset2, m_depth, m_chuteHeight - m_chuteThickness + 0.025f);

	GeometryID obliquePlane = addPlane(1.f, 0.f, 1.f, -m_chuteWidth - m_gapFixOffset, FT_UNFILL);
	setIntersectionType(obliquePlane, IT_SUBTRACT);

	GeometryID baseBoundary = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth - 0.1f, 0.f, -0.125f),
		m_waterBoxWidth - m_chuteWidth + 0.2f, m_depth, 0.1f);

	GeometryID leftBoundary = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(-0.125f, 0.f, m_chuteHeight - 0.1f),
		0.1f, m_depth, 0.2f);
	
	GeometryID rightBoundary = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_waterBoxWidth + 0.025f, 0.f, 0.f),
		0.1f, m_depth, m_waterBoxHeight);

	// boolean value to fill the chute with with water
	const int fillWater = get_option("fillWater", true);

	if(fillWater)
	{
		GeometryID waterBox = addBox(GT_FLUID, FT_SOLID, Point(0.f, 0.f, 0.f),
			m_waterBoxWidth, m_depth, m_waterBoxHeight);
	}
	
	Point bulkPoint = Point(m_chuteUpperWidth, 0.f, m_chuteHeight - m_bulkHeight);
	GeometryID bulk = addBox(GT_FLUID, FT_SOLID, bulkPoint,
		m_bulkWidth, m_depth, m_bulkHeight);
	
	GeometryID obliquePlane2 = addPlane(1.f, 0.f, 1.f, -m_chuteWidth + m_gapFixOffset, FT_UNFILL);
	setEraseOperation(obliquePlane2, ET_ERASE_FLUID);
	setIntersectionType(obliquePlane2, IT_INTERSECT);

	if(fillWater)
	{
		GeometryID fixWaterBox = addBox(GT_FLUID, FT_SOLID, Point(0.f, 0.f, m_chuteHeight + 0.025),
			m_waterBoxWidth, m_depth, m_waterBoxHeight - m_chuteHeight - m_gapFixOffset);
	}

	setPositioning(PP_CENTER);

	// boolean value to a "Tree" (cylinder) in the middle of the chute
	const int addTree = get_option("tree", false);
	
	if(addTree){
		GeometryID tree = addCylinder(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth, m_depth/2.f, m_waterBoxHeight/2.f), 0.1f, m_waterBoxHeight);
	}

	// boolean value to add three "Pillars" (cylinder) at the end of the chute
	const int addPillars = get_option("pillars", false);

	if(addPillars)
	{
		GeometryID leftPillar = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth + 1.f, m_depth/3.f, m_waterBoxHeight/2.f), 0.1f, 0.1f, m_waterBoxHeight);
		GeometryID rightPillar = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth + 1.f, 2*m_depth/3.f, m_waterBoxHeight/2.f), 0.1f, 0.1f, m_waterBoxHeight);
	}

}
