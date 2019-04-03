#include <string>

#include "SubmarineLandslide.h"
#include "GlobalData.h"
#include "cudasimframework.cu"

SubmarineLandslide::SubmarineLandslide(GlobalData *_gdata) :
	XProblem(_gdata),
	// boolean value to fill the chute with with water
	m_fillWater(get_option("m_fillWater", true)),
	// boolean value to a "Tree" (cylinder) in the middle of the chute
	m_addTree(get_option("tree", false)),
	// boolean value to add three "Pillars" (cylinder) at the end of the chute
	m_addPillars(get_option("pillars", false))
{
	m_name = "SubmarineLandslide";

	setFramework();

	const int mlsIters = get_option("mls", 0); // --mls N to enable MLS filter every N iterations
	const int ppH = get_option("ppH", 32); // --ppH N to change deltap to H/N

	if (mlsIters > 0) addFilter(MLS_FILTER, mlsIters);

	set_deltap(m_bulkWidth/ppH);

	m_depth = round_up(m_depth, m_deltap);

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
		periodicity<PERIODIC_Y>,
		add_flags<ENABLE_MULTIFLUID>
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

	m_sand = add_fluid(m_sandDensity);
	set_equation_of_state(m_sand, 7.0f, c0);
	set_dynamic_visc(m_sand, 1.f);
	set_yield_strength(m_sand, m_sandYieldStrength);

	m_water = add_fluid(m_waterDensity);
	set_equation_of_state(m_water, 7.0f, c0);
	set_dynamic_visc(m_water, 1.e-3f);
	set_yield_strength(m_water, m_waterYieldStrength);

	physparams()->dcoeff = 5.0f*g*1;
	physparams()->r0 = m_deltap;
}

bool SubmarineLandslide::is_sand(double4 particlePosition)
{
	double bulkX = particlePosition.x - m_chuteUpperWidth;
	double bulkZ = particlePosition.z; //- ( m_chuteHeight - m_bulkHeight );
	return ( bulkX > 0  && bulkX <= m_bulkWidth && bulkZ <= m_chuteHeight );
}

void SubmarineLandslide::initializeParticles(BufferList &buffers, const uint numParticles)
{
	// 1. warn the user if this is expected to take much time
	printf("Initializing particles density and mass...\n");

	// 2. grab the particle arrays from the buffer list
	float4 *vel = buffers.getData<BUFFER_VEL>();
	particleinfo *info = buffers.getData<BUFFER_INFO>();
	double4 *pos_global = buffers.getData<BUFFER_POS_GLOBAL>();
	float4 *pos = buffers.getData<BUFFER_POS>();

	// 3. iterate on the particles
	const float z_intf = m_chuteHeight;
	const float z_freeSurface = m_fillWater ? m_waterBoxHeight : z_intf;
	// pressure at interface, from heavy fluid
	const float g = length(physparams()->gravity);

	for (uint i = 0; i < numParticles; i++) {
		float rho = 1;
		double depth = z_freeSurface - pos_global[i].z;
		// for boundary particles, we use the density of sand,
		// fluid particles will override fluid_idx depending on whether they are water or sand
		int fluid_idx = m_sand;
		if (FLUID(info[i])) {
			fluid_idx = is_sand(pos_global[i]) ? m_sand : m_water;
			// hydrostatic density: for the heavy fluid, this is simply computed
			// as the density that gives pressure rho g h, with h depth
			rho = hydrostatic_density(depth, fluid_idx);
			// more complex way:
			if (fluid_idx == m_sand) {
				float P = physparams()->rho0[m_water]*(m_waterBoxHeight-z_intf)*g;
				// plus hydrostatic pressure from _our_ fluid
				P += physparams()->rho0[m_sand]*(z_intf - pos_global[i].z)*g;
				rho = density_for_pressure(P, m_sand);
			}
			info[i]= make_particleinfo(PT_FLUID, fluid_idx, i);
		} else if (BOUNDARY(info[i])) {
			rho = hydrostatic_density(depth, fluid_idx);
			info[i]= make_particleinfo(PT_BOUNDARY, fluid_idx, i);
		}
		// fix up the particle mass according to the actual density
		pos[i].w *= physical_density(rho, fluid_idx);
		vel[i].w = rho;
	}	
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

	if(m_fillWater)
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

	if(m_fillWater)
	{
		GeometryID fixWaterBox = addBox(GT_FLUID, FT_SOLID, Point(0.f, 0.f, m_chuteHeight + 0.025),
			m_waterBoxWidth, m_depth, m_waterBoxHeight - m_chuteHeight - m_gapFixOffset);
	}

	setPositioning(PP_CENTER);
	
	if(m_addTree){
		GeometryID tree = addCylinder(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth, m_depth/2.f, m_waterBoxHeight/2.f), 0.1f, m_waterBoxHeight);
	}

	if(m_addPillars)
	{
		GeometryID leftPillar = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth + 1.f, m_depth/3.f, m_waterBoxHeight/2.f), 0.1f, 0.1f, m_waterBoxHeight);
		GeometryID rightPillar = addBox(GT_FIXED_BOUNDARY, FT_SOLID, Point(m_chuteWidth + 1.f, 2*m_depth/3.f, m_waterBoxHeight/2.f), 0.1f, 0.1f, m_waterBoxHeight);
	}

}
