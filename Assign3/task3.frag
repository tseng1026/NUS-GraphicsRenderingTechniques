//============================================================
// STUDENT NAME: TSENG, YU-TING
// STUDENT NO.: A0212195L
// NUS EMAIL ADDRESS: E0503474@u.nus.edu
// COMMENTS TO GRADER:
//
// for CONE, let v  be Ro + t * Rd,
//               Cx be the axis of the cone,
//               with tips of the cone at the origin.
//
// dot(v, Cx)    = dis(v)    * dis(Cx)     * costheta
// dot(v, Cx)**2 = dot(v, v) * dot(Cx, Cx) * costheta**2
//
// dot(v, cx)**2 = dot(Ro + t * Rd, axis)**2
//               = dot(Ro, Cx)**2  + 2t * dot(Ro, Cx)*dot(Rd, Cx) + t**2 * dot(Rd, Cx)**2
//
// dot(v, v)     = dot(Ro + t * Rd, Ro + t * Rd)
//               = dot(Ro, Ro)     + 2t * dot(Ro, Rd)             + t**2 * dot(Rd, Rd)
//
// define cosa equals to dot(Cx, Cx) * costheta
// transform into at**2 + bt + c = 0 format, we get
// a =  dot(Rd, Cx) * dot(Rd, Cx) - dot(Rd, Rd) * cosa**2
// b = (dot(Rd, Cx) * dot(Ro, Cx) - dot(Rd, Ro) * cosa**2) * 2.0
// c =  dot(Ro, Cx) * dot(Ro, Cx) - dot(Ro, Ro) * cosa**2
//
// compute the intersection of ray and cone
//
// eliminate the situation that the dis(v, tips) is larger than the threshold h.
// ============================================================

// FRAGMENT SHADER FOR SHADERTOY
// Run this at https://www.shadertoy.com/new
// See documentation at https://www.shadertoy.com/howto

// Your browser must support WebGL 2.0.
// Check your browser at https://webglreport.com/?v=2


//============================================================================
// Constants.
//============================================================================

const float PI = 3.1415926536;

const vec3 BACKGROUND_COLOR = vec3( 0.1, 0.2, 0.6 );

// Vertical field-of-view angle of camera. In degrees.
const float FOVY = 50.0;

// Use this for avoiding the "epsilon problem" or the shadow acne problem.
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;

// Equivalent to number of recursion levels (0 means ray-casting only).
// We are using iterations to replace recursions.
const int NUM_ITERATIONS = 2;

// Constants for the scene objects.
const int NUM_LIGHTS = 2;
const int NUM_MATERIALS = 4;
const int NUM_PLANES = 2;
const int NUM_SPHERES = 4;
const int NUM_CONES = 2;


//============================================================================
// Define new struct types.
//============================================================================
struct Ray_t {
    vec3 o;  // Ray Origin.
    vec3 d;  // Ray Direction. A unit vector.
};

struct Plane_t {
    // The plane equation is Ax + By + Cz + D = 0.
    float A, B, C, D;
    int materialID;
};

struct Sphere_t {
    vec3 center;
    float radius;
    int materialID;
};

struct Cone_t {
    vec3 tips;      // Tip position.
    vec3 axis;      // Axis.
    float cosa;     // Angle * len(axis).
    float height;
    int materialID;
};

struct Light_t {
    vec3 position;  // Point light 3D position.
    vec3 I_a;       // For Ambient.
    vec3 I_source;  // For Diffuse and Specular.
    vec3 direction;
};

struct Material_t {
    vec3 k_a;   // Ambient coefficient.
    vec3 k_d;   // Diffuse coefficient.
    vec3 k_r;   // Reflected specular coefficient.
    vec3 k_rg;  // Global reflection coefficient.
    float n;    // The specular reflection exponent. Ranges from 0.0 to 128.0.
};

//----------------------------------------------------------------------------
// The lighting model used here is similar to that on Slides 8 and 12 of
// Lecture Topic 9 (Ray Tracing). Here it is computed as
//
//     I_local = SUM_OVER_ALL_LIGHTS {
//                   I_a * k_a +
//                   k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ]
//               }
// and
//     I = I_local  +  k_rg * I_reflected
//----------------------------------------------------------------------------


//============================================================================
// Global scene data.
//============================================================================
Plane_t Plane[NUM_PLANES];
Sphere_t Sphere[NUM_SPHERES];
Cone_t Cone[NUM_CONES];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];



/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////
void InitScene()
{
    // Horizontal plane.
    Plane[0].A = 0.0;
    Plane[0].B = 1.0;
    Plane[0].C = 0.0;
    Plane[0].D = 0.0;
    Plane[0].materialID = 0;

    // Vertical plane.
    Plane[1].A = 0.0;
    Plane[1].B = 0.0;
    Plane[1].C = 1.0;
    Plane[1].D = 3.5;
    Plane[1].materialID = 0;

    // Condition for sphere.
    float offset = 0.5 / sqrt(2.0);
    float current = cos(iTime);
    
    // Right sphere.
    if (current >= 0.0) Sphere[0].center = vec3( -offset * 3.00 - current, 0.8, offset * 3.0 + current );
    if (current <= 0.0) Sphere[0].center = vec3( -offset * 3.00 ,          0.8, offset * 3.0           );
    Sphere[0].radius = 0.5;
    Sphere[0].materialID = 1;
    
    // Left sphere.
    if (current >= 0.0) Sphere[1].center = vec3( offset * 3.00 ,          0.8, -offset * 3.0           );
    if (current <= 0.0) Sphere[1].center = vec3( offset * 3.00 - current, 0.8, -offset * 3.0 + current );
    Sphere[1].radius = 0.5;
    Sphere[1].materialID = 2;
    
    // Middle sphere1.
    Sphere[2].center = vec3( -offset, 0.8, offset );
    Sphere[2].radius = 0.5;
    Sphere[2].materialID = 0;
    
    // Middle sphere2.
    Sphere[3].center = vec3( offset, 0.8, -offset );
    Sphere[3].radius = 0.5;
    Sphere[3].materialID = 0;
    
    // Lamp cone1.
    Cone[0].tips = vec3(-1.0, 3.5, 1.2);
    Cone[0].axis = vec3(0.5 * cos(0.5 * iTime), -1.0, 0.5 * sin(0.5 * iTime));
    Cone[0].cosa = 1.0;
    Cone[0].height = 1.0;
    Cone[0].materialID = 3;
    
    // Lamp cone2.
    Cone[1].tips = vec3(1.2, 3.5, -1.0);
    Cone[1].axis = vec3(0.5 * sin(0.5 * iTime), -1.0, 0.5 * cos(0.5 * iTime));
    Cone[1].cosa = 1.0;
    Cone[1].height = 1.0;
    Cone[1].materialID = 3;

    // Silver material.
    Material[0].k_d = vec3( 0.5, 0.5, 0.5 );
    Material[0].k_a = 0.2 * Material[0].k_d;
    Material[0].k_r = 2.0 * Material[0].k_d;
    Material[0].k_rg = 0.5 * Material[0].k_r;
    Material[0].n = 64.0;

    // Gold material.
    Material[1].k_d = vec3( 0.8, 0.7, 0.1 );
    Material[1].k_a = 0.2 * Material[1].k_d;
    Material[1].k_r = 2.0 * Material[1].k_d;
    Material[1].k_rg = 0.5 * Material[1].k_r;
    Material[1].n = 64.0;

    // Green plastic material.
    Material[2].k_d = vec3( 0.0, 0.8, 0.0 );
    Material[2].k_a = 0.2 * Material[2].k_d;
    Material[2].k_r = vec3( 1.0, 1.0, 1.0 );
    Material[2].k_rg = 0.5 * Material[2].k_r;
    Material[2].n = 128.0;
    
    // Dark blue material
    Material[3].k_d = vec3( 0.1, 0.1, 0.3 );
    Material[3].k_a = 0.5 * Material[3].k_d;
    Material[3].k_r = 2.0 * Material[3].k_d;
    Material[3].k_rg = 0.5 * Material[3].k_r;
    Material[3].n = 128.0;

    // Light 0.
    Light[0].position = Cone[0].tips + Cone[0].axis + vec3(0.0, 0.5, 0.0);
    Light[0].I_a = vec3( 0.1, 0.1, 0.1 );
    Light[0].I_source = vec3( 1.0, 1.0, 1.0 );
    Light[0].direction = vec3( 0.0, 0.0, 0.0 );

    // Light 1.
    Light[1].position = Cone[1].tips + Cone[1].axis + vec3(0.0, 0.5, 0.0);
    Light[1].I_a = vec3( 0.1, 0.1, 0.1 );
    Light[1].I_source = vec3( 1.0, 1.0, 1.0 );
    Light[1].direction = vec3( 0.0, 0.0, 0.0 );
}


/////////////////////////////////////////////////////////////////////////////
// Returns a random number between 0 and 1.
//
// This pseudorandom number generator is based on the 32-bit combined LFSR
// generator proposed in the paper "Tables of Maximally-Equidistributed
// Combined LFSR Generators" by Pierre L'Ecuyer.
// (http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.3639)
/////////////////////////////////////////////////////////////////////////////

// VERY IMPORTANT: The initial seeds rand_z1, rand_z2, rand_z3, rand_z4
// must be larger than 1, 7, 15, and 127 respectively.
const uint CONST_RAND_SEED = 987654321U;
uint rand_z1 = uint(CONST_RAND_SEED + 2U);
uint rand_z2 = uint(CONST_RAND_SEED + 8U);
uint rand_z3 = uint(CONST_RAND_SEED + 16U);
uint rand_z4 = uint(CONST_RAND_SEED + 128U);

float rand(void)
{
    uint b  = ((rand_z1 << 6) ^ rand_z1) >> 13;
    rand_z1 = ((rand_z1 & 4294967294U) << 18) ^ b;
    b       = ((rand_z2 << 2) ^ rand_z2) >> 27;
    rand_z2 = ((rand_z2 & 4294967288U) << 2) ^ b;
    b       = ((rand_z3 << 13) ^ rand_z3) >> 21;
    rand_z3 = ((rand_z3 & 4294967280U) << 7) ^ b;
    b       = ((rand_z4 << 3) ^ rand_z4) >> 12;
    rand_z4 = ((rand_z4 & 4294967168U) << 13) ^ b;
    return float(rand_z1 ^ rand_z2 ^ rand_z3 ^ rand_z4) * 2.3283064365386963e-10;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is such an intersection, outputs the value of t, the position
// of the intersection (hitPos) and the normal vector at the intersection
// (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray,
                     in float tmin, in float tmax,
                     out float t, out vec3 hitPos, out vec3 hitNormal )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;

    // We have a hit -- output results.
    t = t0;
    hitPos = ray.o + t0 * ray.d;
    hitNormal = normalize( N );
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray,
                     in float tmin, in float tmax )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the
// smaller t, the position of the intersection (hitPos) and the normal
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray,
                      in float tmin, in float tmax,
                      out float t, out vec3 hitPos, out vec3 hitNormal )
{
    ///////////////////////////////////
    // TASK 1: WRITE YOUR CODE HERE. //
    // find the coefficient of the quadratic equation
    vec3 Rd = ray.d;
    vec3 Ro = ray.o - sph.center;
    float a = dot(Rd, Rd);         // dot(Rd, Rd) is actually 1
    float b = dot(Rd, Ro) * 2.0;
    float c = dot(Ro, Ro) - sph.radius * sph.radius;
    float D = b * b - 4.0 * a * c;

    // use D to check if there is solution
    // use b to check if the ray direction is toward the sphere
    if (D >= 0.0 && b <= 0.0){
        t = (-b - sqrt(D)) / (2.0 * a);
        hitPos = ray.o + t * ray.d;
        hitNormal = normalize(hitPos - sph.center);
        return (t >= tmin && t <= tmax);
    }
    ///////////////////////////////////

    return false;  // Replace this with your code.

}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray,
                      in float tmin, in float tmax )
{
    ///////////////////////////////////
    // TASK 1: WRITE YOUR CODE HERE. //
    // find the coefficient of the quadratic equation
    vec3 Rd = ray.d;
    vec3 Ro = ray.o - sph.center;
    float a = dot(Rd, Rd);         // dot(Rd, Rd) is actually 1
    float b = dot(Rd, Ro) * 2.0;
    float c = dot(Ro, Ro) - sph.radius * sph.radius;
    float D = b * b - 4.0 * a * c;

    // use D to check if there is solution
    // use b to check if the ray direction is toward the sphere
    if (D >= 0.0 && b <= 0.0){
        float t = (-b - sqrt(D)) / (2.0 * a);
        return (t >= tmin && t<= tmax);
    }
    ///////////////////////////////////

    return false;  // Replace this with your code.

}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a cone and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the
// smaller t, the position of the intersection (hitPos) and the normal
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectCone( in Cone_t cone, in Ray_t ray, in float tmin, in float tmax,
                     out float t, out vec3 hitPos, out vec3 hitNormal )
{
    vec3 Cx = cone.axis;
    vec3 Rd = ray.d;
    vec3 Ro = ray.o - cone.tips;
    float a =  dot(Rd, Cx) * dot(Rd, Cx) - dot(Rd, Rd) * cone.cosa * cone.cosa;  // dot(Rd, Rd) is actually 1
    float b = (dot(Rd, Cx) * dot(Ro, Cx) - dot(Rd, Ro) * cone.cosa * cone.cosa) * 2.0;
    float c =  dot(Ro, Cx) * dot(Ro, Cx) - dot(Ro, Ro) * cone.cosa * cone.cosa;
    float D = b * b - 4.0 * a * c;
    
    // use D to check if there is solution
    if (D >= 0.0){
        float t1 = (-b - sqrt(D)) / (2.0 * a);
        float t2 = (-b + sqrt(D)) / (2.0 * a);
        
        if      (t1 >= tmin && t1 <= tmax) t = t1;
        else if (t2 >= tmin && t2 <= tmax) t = t2;
        else return false;
        
        // relative position
        vec3 tmpPos = Ro + t * Rd;
        
        // only consider one side of the cone, and skip the point out of boundary
        if (dot(tmpPos, Cx) < 0.0) return false;
        if (dot(tmpPos, tmpPos) * cone.cosa > cone.height) return false;
        float height = dot(Cx, tmpPos);
        
        hitPos = ray.o + t * ray.d;
        hitNormal = normalize(tmpPos * height / dot(tmpPos, tmpPos) - cone.axis);
        return (t >= tmin && t<= tmax);
    }

    return false;
    
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a cone and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectCone( in Cone_t cone, in Ray_t ray, in float tmin, in float tmax )
{
    vec3 Cx = cone.axis;
    vec3 Rd = ray.d;
    vec3 Ro = ray.o - cone.tips;
    float a =  dot(Rd, Cx) * dot(Rd, Cx) - dot(Rd, Rd) * cone.cosa * cone.cosa;  // dot(Rd, Rd) is actually 1
    float b = (dot(Rd, Cx) * dot(Ro, Cx) - dot(Rd, Ro) * cone.cosa * cone.cosa) * 2.0;
    float c =  dot(Ro, Cx) * dot(Ro, Cx) - dot(Ro, Ro) * cone.cosa * cone.cosa;
    float D = b * b - 4.0 * a * c;
    
    // use D to check if there is solution
    if (D >= 0.0){
        float t = 0.0;
        float t1 = (-b - sqrt(D)) / (2.0 * a);
        float t2 = (-b + sqrt(D)) / (2.0 * a);
       
        if      (t1 >= tmin && t1 <= tmax) t = t1;
        else if (t2 >= tmin && t2 <= tmax) t = t2;
        else return false;
        
        // relative position
        vec3 tmpPos = Ro + t * Rd;
        
        // only consider one side of the cone, and skip the point out of boundary
        if (dot(tmpPos, Cx) < 0.0) return false;
        if (dot(tmpPos, tmpPos) * cone.cosa > cone.height) return false;
        
        return (t >= tmin && t<= tmax);
    }

    return false;
    
}



/////////////////////////////////////////////////////////////////////////////
// Computes (I_a * k_a) + k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ].
// Input vectors L, N and V are pointing AWAY from surface point.
// Assume all vectors L, N and V are unit vectors.
/////////////////////////////////////////////////////////////////////////////
vec3 PhongLighting( in vec3 L, in vec3 N, in vec3 V, in bool inShadow,
                    in Material_t mat, in Light_t light )
{
    if ( inShadow ) {
        return light.I_a * mat.k_a;
    }
    else {
        vec3 R = reflect( -L, N );
        float N_dot_L = max( 0.0, dot( N, L ) );
        float R_dot_V = max( 0.0, dot( R, V ) );
        float R_dot_V_pow_n = ( R_dot_V == 0.0 )? 0.0 : pow( R_dot_V, mat.n );

        return light.I_a * mat.k_a +
               light.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
    }
}



/////////////////////////////////////////////////////////////////////////////
// Casts a ray into the scene and returns color computed at the nearest
// intersection point. The color is the sum of light from all light sources,
// each computed using Phong Lighting Model, with consideration of
// whether the interesection point is being shadowed from the light.
// If there is no interesection, returns the background color, and outputs
// hasHit as false.
// If there is intersection, returns the computed color, and outputs
// hasHit as true, the 3D position of the intersection (hitPos), the
// normal vector at the intersection (hitNormal), and the k_rg value
// of the material of the intersected object.
/////////////////////////////////////////////////////////////////////////////
vec3 CastRay( in Ray_t ray,
              out bool hasHit, out vec3 hitPos,
              out vec3 hitNormal, out vec3 k_rg )
{
    // Find whether and where the ray hits some object.
    // Take the nearest hit point.

    bool hasHitSomething = false;
    float nearest_t = DEFAULT_TMAX;   // The ray parameter t at the nearest hit point.
    vec3 nearest_hitPos;              // 3D position of the nearest hit point.
    vec3 nearest_hitNormal;           // Normal vector at the nearest hit point.
    int nearest_hitMatID;             // MaterialID of the object at the nearest hit point.

    float temp_t;
    vec3 temp_hitPos;
    vec3 temp_hitNormal;
    bool temp_hasHit;

    ///////////////////////////////////////////////////////////////////////////
    // TASK 1:
    // * Try interesecting input ray with all the planes and spheres,
    //   and record the front-most (nearest) interesection.
    // * If there is interesection, need to record hasHitSomething,
    //   nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////
    // TASK 1: WRITE YOUR CODE HERE. //
    // Try interesecting input ray with all the planes,
    // and record the front-most (nearest) interesection.
    for (int k = 0; k < NUM_PLANES; k ++){
        temp_hasHit = IntersectPlane(Plane[k], ray, DEFAULT_TMIN, DEFAULT_TMAX,
                      temp_t, temp_hitPos, temp_hitNormal);
        
        // If there is interesection, need to record hasHitSomething,
        // nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
        if (temp_hasHit && temp_t < nearest_t){
            hasHitSomething = true;
            nearest_t = temp_t;
            nearest_hitPos = temp_hitPos;
            nearest_hitNormal = temp_hitNormal;
            nearest_hitMatID  = Plane[k].materialID;
        }
    }
    
    // Try interesecting input ray with all the spheres,
    // and record the front-most (nearest) interesection.
    for (int k = 0; k < NUM_SPHERES; k ++){
        temp_hasHit = IntersectSphere(Sphere[k], ray, DEFAULT_TMIN, DEFAULT_TMAX,
                      temp_t, temp_hitPos, temp_hitNormal);
        
        // If there is interesection, need to record hasHitSomething,
        // nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
        if (temp_hasHit && temp_t < nearest_t){
            hasHitSomething = true;
            nearest_t = temp_t;
            nearest_hitPos = temp_hitPos;
            nearest_hitNormal = temp_hitNormal;
            nearest_hitMatID  = Sphere[k].materialID;
        }
    }
    
    // Try interesecting input ray with all the cones,
    // and record the front-most (nearest) interesection.
    for (int k = 0; k < NUM_CONES; k ++){
        temp_hasHit = IntersectCone(Cone[k], ray, DEFAULT_TMIN, DEFAULT_TMAX,
                      temp_t, temp_hitPos, temp_hitNormal);
        
        // If there is interesection, need to record hasHitSomething,
        // nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
        if (temp_hasHit && temp_t < nearest_t){
            hasHitSomething = true;
            nearest_t = temp_t;
            nearest_hitPos = temp_hitPos;
            nearest_hitNormal = temp_hitNormal;
            nearest_hitMatID  = Cone[k].materialID;
        }
    }
    ///////////////////////////////////


    // One of the output results.
    hasHit = hasHitSomething;
    if ( !hasHitSomething ) return BACKGROUND_COLOR;

    vec3 I_local = vec3( 0.0 );  // Result color will be accumulated here.

    ///////////////////////////////////////////////////////////////////////////
    // TASK 1:
    // * Accumulate lighting from each light source on the nearest hit point.
    //   They are all accumulated into I_local.
    // * For each light source, make a shadow ray, and check if the shadow ray
    //   intersects any of the objects (the planes and spheres) between the
    //   nearest hit point and the light source.
    // * Then, call PhongLighting() to compute lighting for this light source.
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////
    // TASK 1: WRITE YOUR CODE HERE. //
    // Accumulate lighting from each light source on the nearest hit point.
    Ray_t ShadowRay[NUM_LIGHTS];
    bool inshadow[NUM_LIGHTS];
    for (int l = 0; l < NUM_LIGHTS; l ++){
        // For each light source, make a shadow ray
        ShadowRay[l].o = nearest_hitPos;
        ShadowRay[l].d = normalize(Light[l].position - ShadowRay[l].o);

        // check if the shadow ray intersects any of the planes
        // between the nearest hit point and the light source.
        for (int k = 0; k < NUM_PLANES; k ++){
            if(inshadow[l]) break;
            inshadow[l] = IntersectPlane(Plane[k], ShadowRay[l], DEFAULT_TMIN,
                          distance(nearest_hitPos, Light[l].position));
        }
        
        // check if the shadow ray intersects any of the spheres
        // between the nearest hit point and the light source.
        for (int k = 0; k < NUM_SPHERES; k ++){
            if(inshadow[l]) break;
            inshadow[l] = IntersectSphere(Sphere[k], ShadowRay[l], DEFAULT_TMIN,
                          distance(nearest_hitPos, Light[l].position));
        }
        
        // check if the shadow ray intersects any of the cones
        // between the nearest hit point and the light source.
        for (int k = 0; k < NUM_CONES; k ++){
            if(inshadow[l]) break;
            inshadow[l] = IntersectCone(Cone[k], ShadowRay[l], DEFAULT_TMIN,
                          distance(nearest_hitPos, Light[l].position));
        }

        // Then, call PhongLighting() to compute lighting for this light source.
        I_local += PhongLighting(ShadowRay[l].d, nearest_hitNormal, -ray.d, inshadow[l],
                   Material[nearest_hitMatID], Light[l]);
    }
    ///////////////////////////////////

    // Populate output results.
    hitPos = nearest_hitPos;
    hitNormal = nearest_hitNormal;
    k_rg = Material[nearest_hitMatID].k_rg;

    return I_local;
}



/////////////////////////////////////////////////////////////////////////////
// Execution of fragment shader starts here.
// 1. Initializes the scene.
// 2. Compute a primary ray for the current pixel (fragment).
// 3. Trace ray into the scene with NUM_ITERATIONS recursion levels.
/////////////////////////////////////////////////////////////////////////////
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Initialize random number generator before the first call to rand().
    uint RAND_SEED = uint( (mod(iTime*100.0, 100.0) + 101.01) *
                           (fragCoord.x + 17.0) * (fragCoord.y + 23.0) );
    rand_z1 = uint(RAND_SEED + 2U);
    rand_z2 = uint(RAND_SEED + 8U);
    rand_z3 = uint(RAND_SEED + 16U);
    rand_z4 = uint(RAND_SEED + 128U);
    
    InitScene();


    // Camera position and orientation in world space.
    vec3 cam_pos = vec3( 3.0, 1.0, 3.0 );
    vec3 cam_lookat = vec3( 0.0, 1.2, 0.0);
    vec3 cam_up_vec = vec3( 0.0, 1.0, 0.0 );

    // Camera coordinate frame in world space.
    vec3 cam_z_axis = normalize( cam_pos - cam_lookat );
    vec3 cam_x_axis = normalize( cross(cam_up_vec, cam_z_axis) );
    vec3 cam_y_axis = normalize( cross(cam_z_axis, cam_x_axis));

    // Vertical field-of-view angle of camera. In radians.
    float cam_FOVY = FOVY * PI / 180.0;

    // Perpendicular distance of the image rectangle from the camera.
    // If implementing depth-of-field, the plane of the image rectangle
    // is the plane of focus.
    float image_dist = distance(cam_pos, vec3(0.0, 0.7, 0.0));

    float image_height = 2.0 * image_dist * tan(cam_FOVY / 2.0);
    float image_width = image_height * iResolution.x / iResolution.y;
    float pixel_width = image_width / iResolution.x;

    // Image rectangle origin (bottom-leftmost corner) position in camera space.
    vec3 image_origin = vec3(-image_width/2.0, -image_height/2.0, -image_dist);


    /////////////////////////////////////////////////////////////////////////
    // TASK 2:
    // * Trace multiple (SPP) random primary rays per pixel to produce
    //   depth-of-field effect and for image anti-aliasing (reduce jaggies).
    // * Each primary ray starts from a random position on the lens and
    //   points towards a random position inside the current pixel.
    // * The lens is assumed to have a square-shaped aperture of size
    //   aperture_width x aperture_width. The lens is centered at the
    //   the origin of the camera frame, and parallel to the x-y plane
    //   of the camera frame.
    // * The final color of the current pixel is the average color of
    //   all the primary rays.
    /////////////////////////////////////////////////////////////////////////

    //=======================================================================
    // These constants are used for distribution ray tracing to produce
    // depth-of-field effect and for image anti-aliasing (reduce jaggies).
    //=======================================================================
    // Number of samples (random primary rays) per pixel.
    const int SPP = 32;

    // Lens aperture width. Assume square aperture.
    const float aperture_width = 0.3;
    //=======================================================================

    ////////////////////////////////////
    // TASK 2: MODIFY THE CODE BELOW. //
    vec3 I_result = vec3( 0.0 );

    for (int k = 0; k < SPP; k ++){
        float rand_pixel_x = pixel_width * (rand() - 0.5) + pixel_width * fragCoord.x;
        float rand_pixel_y = pixel_width * (rand() - 0.5) + pixel_width * fragCoord.y;
        
        vec3 rand_pixel_pos = image_origin +
                                  vec3(rand_pixel_x, rand_pixel_y, 0.0);
        
        // Create random primary ray.
        float rand_cam_x = (rand() - 0.5) * aperture_width;
        float rand_cam_y = (rand() - 0.5) * aperture_width;
        
        Ray_t randRay;
        randRay.o = cam_pos + vec3(rand_cam_x, rand_cam_y, 0);
        randRay.d = normalize( cam_pos + rand_pixel_pos.x * cam_x_axis  +
                                         rand_pixel_pos.y * cam_y_axis  +
                                         rand_pixel_pos.z * cam_z_axis - randRay.o);

        // Start Ray Tracing.
        // Use iterations to emulate the recursion.

        vec3 compounded_k_rg = vec3( 1.0 );
        Ray_t nextRay = randRay;

        for ( int level = 0; level <= NUM_ITERATIONS; level++ )
        {
            bool hasHit;
            vec3 hitPos, hitNormal, k_rg;

            vec3 I_local = CastRay( nextRay, hasHit, hitPos, hitNormal, k_rg );

            I_result += compounded_k_rg * I_local;

            if ( !hasHit ) break;

            compounded_k_rg *= k_rg;

            nextRay = Ray_t( hitPos, normalize( reflect(nextRay.d, hitNormal) ) );
        }
    }
    
    ////////////////////////////////////

    fragColor = vec4( I_result / float(SPP), 1.0 );
}
