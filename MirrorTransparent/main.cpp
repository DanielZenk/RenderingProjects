#include "main.hpp"
#include<cmath>

// ******* Function Member Implementation *******

// ***** Shape function members *****
Shape::Shape() : mColour{ 0, 0, 0 }
{}

void Shape::setColour(Colour const& col)
{
    mColour = col;
}

Colour Shape::getColour() const
{
    return mColour;
}

void Shape::setMaterial(std::shared_ptr<Material> const& material)
{
    mMaterial = material;
}

std::shared_ptr<Material> Shape::getMaterial() const
{
    return mMaterial;
}

// ***** Sampler function members *****
Sampler::Sampler(int numSamples, int numSets) :
    mNumSamples{ numSamples }, mNumSets{ numSets }, mCount{ 0 }, mJump{ 0 }
{
    mSamples.reserve(mNumSets* mNumSamples);
    setupShuffledIndeces();
}

int Sampler::getNumSamples() const
{
    return mNumSamples;
}

void Sampler::setupShuffledIndeces()
{
    mShuffledIndeces.reserve(mNumSamples * mNumSets);
    std::vector<int> indices;

    std::random_device d;
    std::mt19937 generator(d());

    for (int j = 0; j < mNumSamples; ++j)
    {
        indices.push_back(j);
    }

    for (int p = 0; p < mNumSets; ++p)
    {
        std::shuffle(indices.begin(), indices.end(), generator);

        for (int j = 0; j < mNumSamples; ++j)
        {
            mShuffledIndeces.push_back(indices[j]);
        }
    }
}

atlas::math::Point Sampler::sampleUnitSquare()
{
    if (mCount % mNumSamples == 0)
    {
        atlas::math::Random<int> engine;
        mJump = ((int)rand() % mNumSets) * mNumSamples;
    }

    return mSamples[mJump + mShuffledIndeces[mJump + mCount++ % mNumSamples]];
}

void Sampler::mapToHemisphere(const float e) {
    size_t size = mSamples.size();
    mHemisphereSamples.reserve(mNumSets * mNumSamples);
    for (int i = 0; i < size; i++) {
        double cosPhi = cos(2.0 * glm::pi<float>() * mSamples[i].x);
        double sinPhi = sin(2.0 * glm::pi<float>() * mSamples[i].x);
        double cosTheta = pow((1.0 * mSamples[i].y), 1.0 / (e + 1.0));
        double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        double pointU = sinTheta * cosPhi;
        double pointV = sinTheta * sinPhi;
        double pointW = cosTheta;
        mHemisphereSamples.push_back(atlas::math::Point(pointU, pointV, pointW));
    }
}

atlas::math::Point Sampler::sampleHemisphere() {
    if (mCount % mNumSamples == 0)
    {
        atlas::math::Random<int> engine;
        mJump = (engine.getRandomMax() % mNumSets) * mNumSamples;
    }

    return mSamples[mJump + mShuffledIndeces[mJump + mCount++ % mNumSamples]];
}

// ***** Light function members *****
Colour Light::L([[maybe_unused]] ShadeRec& sr)
{
    return mRadiance * mColour;
}

void Light::scaleRadiance(float b)
{
    mRadiance = b;
}

void Light::setColour(Colour const& c)
{
    mColour = c;
}

// ***** Sphere function members *****
Sphere::Sphere(atlas::math::Point center, float radius) :
    mCentre{ center }, mRadius{ radius }, mRadiusSqr{ radius * radius }
{}

bool Sphere::hit(atlas::math::Ray<atlas::math::Vector> const& ray,
    ShadeRec& sr) const
{
    atlas::math::Vector tmp = ray.o - mCentre;
    float t{ std::numeric_limits<float>::max() };
    bool intersect{ intersectRay(ray, t) };

    // update ShadeRec info about new closest hit
    if (intersect && t < sr.t)
    {
        sr.normal = (tmp + t * ray.d) / mRadius;
        sr.ray = ray;
        sr.color = mColour;
        sr.t = t;
        sr.material = mMaterial;
    }

    return intersect;
}

bool Sphere::intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray,
    float& tMin) const
{
    const auto tmp{ ray.o - mCentre };
    const auto a{ glm::dot(ray.d, ray.d) };
    const auto b{ 2.0f * glm::dot(ray.d, tmp) };
    const auto c{ glm::dot(tmp, tmp) - mRadiusSqr };
    const auto disc{ (b * b) - (4.0f * a * c) };

    if (atlas::core::geq(disc, 0.0f))
    {
        const float kEpsilon{ 0.01f };
        const float e{ std::sqrt(disc) };
        const float denom{ 2.0f * a };

        // Look at the negative root first
        float t = (-b - e) / denom;
        if (atlas::core::geq(t, kEpsilon))
        {
            tMin = t;
            return true;
        }

        // Now the positive root
        t = (-b + e) / denom;
        if (atlas::core::geq(t, kEpsilon))
        {
            tMin = t;
            return true;
        }
    }

    return false;
}

bool Sphere::shadowHit(atlas::math::Ray<atlas::math::Vector> const& ray, float& t) {
    using atlas::math::Point;
    using atlas::math::Vector;

    Vector temp = ray.o - mCentre;
    float a{ glm::dot(ray.d, ray.d) };
    float b{ glm::dot(ray.d * 2.0f, temp) };
    float c{ glm::dot(temp, temp) - (mRadius * mRadius) };
    float disc{ b * b - 4.0f * a * c };
    if (disc < 0.0f) {
        return false;
    }
    float newT{ (-b - sqrt(disc)) / (2.0f * a) };
    if (newT > 0.1) {
        t = newT;
        return true;
    }
    /*newT = (-b + sqrt(disc)) / (2.0f * a);
    if (newT > 0.1) {
        t = newT;
        return true;
    }*/
    return false;
}

Triangle::Triangle(atlas::math::Point mp0, atlas::math::Point mp1, atlas::math::Point mp2) :
    p1{ mp1 }, p2{ mp2 }, p0{ mp0 }
{}

bool Triangle::hit([[maybe_unused]] atlas::math::Ray<atlas::math::Vector> const& ray, [[maybe_unused]] ShadeRec& sr) const {
    //Get normal of triangle
    using atlas::math::Vector;

    Vector help1 = p1 - p0;
    Vector help2 = p2 - p0;
    Vector triNorm = glm::cross(help1, help2);

    float helper = glm::dot(triNorm, -ray.d);

    float t{ std::numeric_limits<float>::max() };
    bool intersect{ intersectRay(ray, t) };

    // update ShadeRec info about new closest hit
    if (intersect && t < sr.t)
    {
        if (helper > 0) {
            sr.normal = glm::normalize(triNorm);
        }
        else {
            sr.normal = -glm::normalize(triNorm);

        }
        sr.ray = ray;
        sr.color = mColour;
        sr.t = t;
        sr.material = mMaterial;
    }

    return intersect;
}

bool Triangle::intersectRay([[maybe_unused]] atlas::math::Ray<atlas::math::Vector> const& ray,
    [[maybe_unused]] float& tMin) const {
    using atlas::math::Vector;

    Vector help1 = p1 - p0;
    Vector help2 = p2 - p0;
    Vector triNorm = glm::cross(help1, help2);
    float d = glm::dot(p0, triNorm);
    float x = glm::dot(triNorm, ray.d);
    if (x == 0) return false;
    float t = (glm::dot(ray.o, triNorm) + d) / x;
    //if we made it this far, we intersected the triangles plane.
    //next we check to make sure the intersect point is in the triangle
    Vector p = ray.o + t * ray.d;
    Vector edge0 = p1 - p0;
    Vector edge1 = p2 - p1;
    Vector edge2 = p0 - p2;
    Vector c0 = glm::cross(edge0, p - p0);
    Vector c1 = glm::cross(edge1, p - p1);
    Vector c2 = glm::cross(edge2, p - p2);

    if (glm::dot(c0, triNorm) < 0) {
        return false;
    }
    else if (glm::dot(c1, triNorm) < 0) {
        return false;
    }
    else if (glm::dot(c2, triNorm) < 0) {
        return false;
    }
    else {
        tMin = t;
        return true;
    }
}

bool Triangle::shadowHit([[maybe_unused]] atlas::math::Ray<atlas::math::Vector> const& ray, [[maybe_unused]] float& t) {
    using atlas::math::Vector;

    Vector help1 = p1 - p0;
    Vector help2 = p2 - p0;
    Vector triNorm = glm::cross(help1, help2);
    float d = glm::dot(p0, triNorm);
    float x = glm::dot(triNorm, ray.d);
    if (x == 0) return false;
    float t1 = (glm::dot(ray.o, triNorm) + d) / x;
    //if we made it this far, we intersected the triangles plane.
    //next we check to make sure the intersect point is in the triangle
    Vector p = ray.o + t1 * ray.d;
    Vector edge0 = p1 - p0;
    Vector edge1 = p2 - p1;
    Vector edge2 = p0 - p2;
    Vector c0 = glm::cross(edge0, p - p0);
    Vector c1 = glm::cross(edge1, p - p1);
    Vector c2 = glm::cross(edge2, p - p2);

    if (glm::dot(c0, triNorm) < 0) {
        return false;
    }
    else if (glm::dot(c1, triNorm) < 0) {
        return false;
    }
    else if (glm::dot(c2, triNorm) < 0) {
        return false;
    }
    else {
        t = t1;
        return true;
    }
}

Plane::Plane(atlas::math::Point mp0, atlas::math::Vector mnorm) :
    p0{ mp0 }, norm{ mnorm }
{}

bool Plane::hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const {
    using atlas::math::Vector;

    float t{ std::numeric_limits<float>::max() };
    bool intersect{ intersectRay(ray, t) };

    // update ShadeRec info about new closest hit
    if (intersect && t < sr.t)
    {
        sr.normal = -norm;
        sr.ray = ray;
        sr.color = mColour;
        sr.t = t;
        sr.material = mMaterial;
    }

    return intersect;
}

bool Plane::intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray,
    float& tMin) const {

    using atlas::math::Vector;

    Vector temp{ p0 - ray.o };
    float a{ glm::dot(temp, norm) };
    float b{ glm::dot(ray.d, norm) };
    if (b == 0) return false;
    float t{ a / b };
    if (t < 0) return false;
    tMin = t;
    return true;
}

bool Plane::shadowHit([[maybe_unused]] atlas::math::Ray<atlas::math::Vector> const& ray,
    [[maybe_unused]] float& t) {
    return false;
}

// ***** Camera function members *****
Camera::Camera() :
    mEye{ 0.0f, 0.0f, 500.0f },
    mLookAt{ 0.0f },
    mUp{ 0.0f, 1.0f, 0.0f },
    mU{ 1.0f, 0.0f, 0.0f },
    mV{ 0.0f, 1.0f, 0.0f },
    mW{ 0.0f, 0.0f, 1.0f }
{}

void Camera::setEye(atlas::math::Point const& eye)
{
    mEye = eye;
}

void Camera::setLookAt(atlas::math::Point const& lookAt)
{
    mLookAt = lookAt;
}

void Camera::setUpVector(atlas::math::Vector const& up)
{
    mUp = up;
}

void Camera::computeUVW()
{
    mW = glm::normalize(mEye - mLookAt);
    mU = glm::normalize(glm::cross(mUp, mW));
    mV = glm::cross(mW, mU);

    if (areEqual(mEye.x, mLookAt.x) && areEqual(mEye.z, mLookAt.z) &&
        mEye.y > mLookAt.y)
    {
        mU = { 0.0f, 0.0f, 1.0f };
        mV = { 1.0f, 0.0f, 0.0f };
        mW = { 0.0f, 1.0f, 0.0f };
    }

    if (areEqual(mEye.x, mLookAt.x) && areEqual(mEye.z, mLookAt.z) &&
        mEye.y < mLookAt.y)
    {
        mU = { 1.0f, 0.0f, 0.0f };
        mV = { 0.0f, 0.0f, 1.0f };
        mW = { 0.0f, -1.0f, 0.0f };
    }
}

// ***** Multijittered function members *****
Jittered::Jittered(int numSamples, int numSets) : Sampler{ numSamples, numSets }
{
    generateSamples();
}

void Jittered::generateSamples()
{
    int n = static_cast<int>(glm::sqrt(static_cast<float>(mNumSamples)));

    atlas::math::Random<float> engine;

    for (int p = 0; p < mNumSets; p++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                mSamples.push_back(atlas::math::Point{ (k + engine.getRandomOne()) / n, (j + engine.getRandomOne()) / n, 0.0f });
            }
        }
    }

    mapToHemisphere(1.0f);
}

// ***** Pinhole function members *****
Pinhole::Pinhole() : Camera{}, mDistance{ 500.0f }, mZoom{ 1.0f }
{}

void Pinhole::setDistance(float distance)
{
    mDistance = distance;
}

void Pinhole::setZoom(float zoom)
{
    mZoom = zoom;
}

atlas::math::Vector Pinhole::rayDirection(atlas::math::Point const& p) const
{
    const auto dir = p.x * mU + p.y * mV - mDistance * mW;
    return glm::normalize(dir);
}

void Pinhole::renderScene(std::shared_ptr<World> world) const
{
    using atlas::math::Point;
    using atlas::math::Ray;
    using atlas::math::Vector;

    Point samplePoint{}, pixelPoint{};
    Ray<atlas::math::Vector> ray{};

    ray.o = mEye;
    float avg{ 1.0f / world->sampler->getNumSamples() };

    for (int r{ 0 }; r < world->height; ++r)
    {
        for (int c{ 0 }; c < world->width; ++c)
        {
            Colour pixelAverage{ 0, 0, 0 };

            for (int j = 0; j < world->sampler->getNumSamples(); ++j)
            {
                ShadeRec trace_data{};
                trace_data.world = world;
                trace_data.t = std::numeric_limits<float>::max();
                samplePoint = world->sampler->sampleUnitSquare();
                pixelPoint.x = c - 0.5f * world->width + samplePoint.x;
                pixelPoint.y = r - 0.5f * world->height + samplePoint.y;
                ray.d = rayDirection(pixelPoint);
                bool hit{};

                for (auto obj : world->scene)
                {
                    hit |= obj->hit(ray, trace_data);
                }

                if (hit)
                {
                    trace_data.hitPoint = ray.o + trace_data.t * ray.d;
                    pixelAverage += world->tracer->trace_ray(ray, 0, world); //trace_data.material->shade(trace_data);//+ 
                }
            }
            pixelAverage.r *= avg;
            pixelAverage.g *= avg;
            pixelAverage.b *= avg;

            world->image.push_back(outOfGamut(pixelAverage));
        }
    }
}

// ***** Lambertian function members *****
Lambertian::Lambertian() : mDiffuseColour{}, mDiffuseReflection{}
{}

Lambertian::Lambertian(Colour diffuseColor, float diffuseReflection) :
    mDiffuseColour{ diffuseColor }, mDiffuseReflection{ diffuseReflection }
{}

Colour
Lambertian::fn([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected,
    [[maybe_unused]] atlas::math::Vector const& incoming) const
{
    return mDiffuseColour * mDiffuseReflection * glm::one_over_pi<float>();
}

Colour
Lambertian::rho([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected) const
{
    return mDiffuseColour * mDiffuseReflection;
}

void Lambertian::setDiffuseReflection(float kd)
{
    mDiffuseReflection = kd;
}

void Lambertian::setDiffuseColour(Colour const& colour)
{
    mDiffuseColour = colour;
}

Colour Lambertian::sample_f(const ShadeRec& sr,[[maybe_unused]] const Vector& wo, Vector& wi,
                            float& pdf, std::shared_ptr<World> world) {
    Vector w = sr.normal;
    Vector v = { 0.0034, 1.0, 0.0071 };
    v = glm::cross(v, w);
    v = glm::normalize(v);
    Vector u = glm::cross(v, w);


    Point sp = world->sampler->sampleHemisphere();
    wi = Vector{sp.x * u.x, sp.x * u.y, sp.x * u.z} + 
         Vector{sp.y * v.x, sp.y * v.y, sp.y * v.z} +
         Vector{sp.z * w.x, sp.z * w.y, sp.z * w.z};
    wi = glm::normalize(wi);
    pdf =  (glm::dot(sr.normal, wi)) / glm::pi<float>();

    return (mDiffuseReflection * mDiffuseColour * glm::one_over_pi<float>());

}

// ***** Perfect transmitter function members
PerfectTransmitter::PerfectTransmitter() {}

PerfectTransmitter::PerfectTransmitter(Colour mkt, float mior) :
    kt{ mkt }, ior{ mior }
{}

Colour PerfectTransmitter::fn([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected,
    [[maybe_unused]] atlas::math::Vector const& incoming) const {
    return { 1, 1, 1 };
}

Colour PerfectTransmitter::rho([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected) const
{
    return { 1, 1, 1 };
}

void PerfectTransmitter::setDiffuseReflection(float kd)
{
    ior = kd;
}

void PerfectTransmitter::setDiffuseColour(Colour const& colour)
{
    kt = colour;
}

bool PerfectTransmitter::tir(ShadeRec const& sr) const 
{
    Vector wo = -sr.ray.d;
    float cosThetaI = glm::dot(sr.normal, wo);
    float eta = ior;

    if (cosThetaI < 0.0) {
        eta = 1.0f / eta;
    }

    return (1.0f - (1.0f - (cosThetaI * cosThetaI)) / (eta * eta) < 0.0f);
}

Colour PerfectTransmitter::sample_f(const ShadeRec& sr,
    atlas::math::Vector const& reflected,
    atlas::math::Vector& incoming) const 
{
    atlas::math::Normal n = sr.normal;
    float cosThetaI = glm::dot(n, reflected);
    float eta = ior;

    if (cosThetaI < 0.0) {
        cosThetaI = -cosThetaI;
        n = -n;
        eta = 1.0f / eta;
    }

    float temp = 1.0f - (1.0f - (cosThetaI * cosThetaI)) / (eta * eta);
    float cosTheta2 = sqrt(temp);
    float multiplier = (cosTheta2 - cosThetaI / eta);
    incoming = -reflected / eta - multiplier * n;
    Colour white = { 1, 1, 1 };
    return (kt / (eta * eta) * white / fabs(glm::dot(sr.normal, incoming)));
}

// ***** Phong function members *****
Phong::Phong() :
    Material{},
    mDiffuseBRDF{ std::make_shared<Lambertian>() },
    mAmbientBRDF{ std::make_shared<Lambertian>() },
    mSpecularBRDF{ std::make_shared<GlossySpecular>() }
{}

Phong::Phong(float kd, float ka, Colour color, float exp) :
    Material{},
    mDiffuseBRDF{ std::make_shared<Lambertian>() },
    mAmbientBRDF{ std::make_shared<Lambertian>() },
    mSpecularBRDF{ std::make_shared<GlossySpecular>() }
{
    setDiffuseReflection(kd);
    setAmbientReflection(ka);
    setDiffuseColour(color);
    setSpecularExp(exp);
}

void Phong::setSpecularExp(float exp) {
    mSpecularBRDF->setExp(exp);
}

void Phong::setDiffuseReflection(float k) {
    mDiffuseBRDF->setDiffuseReflection(k);
    mSpecularBRDF->setDiffuseReflection(k);
}

void Phong::setAmbientReflection(float k)
{
    mAmbientBRDF->setDiffuseReflection(k);
}

void Phong::setDiffuseColour(Colour c) {
    mDiffuseBRDF->setDiffuseColour(c);
    mAmbientBRDF->setDiffuseColour(c);
    mSpecularBRDF->setDiffuseColour({ 1,1,1 });
}

Colour Phong::shade(ShadeRec& sr) {
    using atlas::math::Ray;
    using atlas::math::Vector;

    Vector wo = glm::normalize(-sr.ray.d);
    Colour L = mAmbientBRDF->rho(sr, wo) * sr.world->ambient->L(sr);
    size_t numLights = sr.world->lights.size();

    for (size_t i{ 0 }; i < numLights; ++i)
    {
        Vector wi = sr.world->lights[i]->getDirection(sr);
        float nDotWi = glm::dot(sr.normal, wi);

        if (nDotWi > 0.0f)
        {

            Ray<atlas::math::Vector> shadowRay{};
            shadowRay.o = sr.hitPoint;
            shadowRay.d = wi;
            bool inShadow = sr.world->lights[i]->inShadow(shadowRay, sr);
            if (!inShadow) {
                L += (mDiffuseBRDF->fn(sr, wo, wi) + mSpecularBRDF->fn(sr, wo, wi)) * sr.world->lights[i]->L(sr) *
                    nDotWi;
            }
        }
    }

    return L;
}

Colour Phong::path_shade([[maybe_unused]] ShadeRec& sr, [[maybe_unused]] std::shared_ptr<World> world) {
    return { 0, 0, 0 };
}

// ***** Reflective function members
Reflective::Reflective(void) : Phong{} {};


Reflective::Reflective(float kd, float ka, Colour color, float exp) :
    Phong{ kd, ka, color, exp },
    mPerfectBRDF{ new PerfectSpecular() } {};

Colour Reflective::shade(ShadeRec& sr) {
    Colour L{ Phong::shade(sr) };

    Vector wo = -sr.ray.d;
    Vector wi;
    Colour fr = mPerfectBRDF->sample_f(sr, wo, wi);

    Ray reflect(sr.hitPoint, wi);
    Colour traceColor = sr.world->tracer->trace_ray(reflect, sr.depth + 1, sr.world);

    Colour returnColor = { fr.r + traceColor.r, fr.g + traceColor.g, fr.b + traceColor.b };
    returnColor *= glm::dot(sr.normal, wi);

    L += returnColor;
    return L;
}

// ***** Reflective function members
Transparent::Transparent(void) : Phong{} {};


Transparent::Transparent(float kd, float ka, Colour color, float exp, Colour kt, float ior) :
    Phong{ kd, ka, color, exp },
    mPerfectBRDF{ new PerfectSpecular() },
    mTransmitterBRDF{ new PerfectTransmitter(kt, ior) }{};

Colour Transparent::shade(ShadeRec& sr) {
    Colour L{ Phong::shade(sr) };

    Vector wo = -sr.ray.d;
    Vector wi;
    Colour fr = mPerfectBRDF->sample_f(sr, wo, wi);

    Ray reflect(sr.hitPoint, wi);

    if (mTransmitterBRDF->tir(sr)) {
        L += sr.world->tracer->trace_ray(reflect, sr.depth + 1, sr.world);
    }
    else {
        Vector wt;
        Colour ft = mTransmitterBRDF->sample_f(sr, wo, wt);
        Ray transmitted(sr.hitPoint, wt);

        Colour traceColor = sr.world->tracer->trace_ray(reflect, sr.depth + 1, sr.world);
        Colour returnColor = fr * traceColor * fabs(glm::dot(sr.normal, wi));

        Colour traceColor2 = sr.world->tracer->trace_ray(transmitted, sr.depth + 1, sr.world);
        Colour returnColor2 = ft * traceColor2 * fabs(glm::dot(sr.normal, wt));

        L += returnColor;
        L += returnColor2;
    }
    return L;
}


Tracer::Tracer() {}

// **** Whitted Tracer function members
WhittedTracer::WhittedTracer()
{
}

Colour WhittedTracer::trace_ray(const Ray ray, const int depth, std::shared_ptr<World> world) const
{
    if (depth > world->maxDepth) {
        return { 0, 0, 0 };
    }
    else {
        ShadeRec sr{};
        sr.t = std::numeric_limits<float>::max();
        bool hit{};
        for (auto obj : world->scene)
        {
            hit |= obj->hit(ray, sr);
        }

        if (hit)
        {
            sr.depth = depth;
            sr.ray = ray;
            sr.hitPoint = ray.o + sr.t * ray.d;
            sr.world = world;
            return sr.material->shade(sr);
        }
        else {
            return world->background;
        }
    }
}

// ***** Glossy Specular function members *****
GlossySpecular::GlossySpecular() :
    mExp{},
    mDiffuseColour{},
    mDiffuseReflection{}
{}

GlossySpecular::GlossySpecular(Colour diffuseColor, float diffuseReflection, float exp) :
    mExp{ exp },
    mDiffuseColour{ diffuseColor },
    mDiffuseReflection{ diffuseReflection }
{}

void GlossySpecular::setExp(float exp) {
    mExp = exp;
}

void GlossySpecular::setDiffuseReflection(float kd)
{
    mDiffuseReflection = kd;
}

void GlossySpecular::setDiffuseColour(Colour const& colour)
{
    mDiffuseColour = colour;
}

Colour GlossySpecular::rho([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected) const
{
    return { 0, 0, 0 };
}

Colour GlossySpecular::fn([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected,
    [[maybe_unused]] atlas::math::Vector const& incoming) const
{
    Colour L = { 0,0,0 };
    float nDotIncoming = glm::dot(sr.normal, incoming);
    atlas::math::Vector r{ -incoming + 2.0f * sr.normal * nDotIncoming };
    //r = glm::normalize(r);
    float rDotOut = glm::dot(r, reflected);
    Colour test = { 1,1,1 };
    if (rDotOut > 0.0) {
        L = test * pow(rDotOut, mExp);
    }
    return L;
}

// ***** Perfect Specular function members *****
PerfectSpecular::PerfectSpecular() :
    mDiffuseColour{},
    mDiffuseReflection{}
{}

PerfectSpecular::PerfectSpecular(Colour diffuseColor, float diffuseReflection) :
    mDiffuseColour{ diffuseColor },
    mDiffuseReflection{ diffuseReflection }
{}

void PerfectSpecular::setDiffuseReflection(float kd)
{
    mDiffuseReflection = kd;
}

void PerfectSpecular::setDiffuseColour(Colour const& colour)
{
    mDiffuseColour = colour;
}

Colour PerfectSpecular::rho([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected) const
{
    return { 0, 0, 0 };
}

Colour PerfectSpecular::fn([[maybe_unused]] ShadeRec const& sr,
    [[maybe_unused]] atlas::math::Vector const& reflected,
    [[maybe_unused]] atlas::math::Vector const& incoming) const
{
    return { 0,0,0 };
}

Colour PerfectSpecular::sample_f(ShadeRec const& sr,
    atlas::math::Vector const& reflected,
    atlas::math::Vector& incoming) const
{
    float nDotWo = glm::dot(sr.normal, reflected);
    atlas::math::Vector test = 2.0f * sr.normal;
    test = test * nDotWo;
    incoming = Vector{ -reflected.x + test.x, -reflected.y + test.y, -reflected.z + test.z };

    return ((mDiffuseReflection * mDiffuseColour) / (glm::dot(sr.normal, incoming)));
    //incoming = -reflected + test;
}

// ***** Matte function members *****
Matte::Matte() :
    Material{},
    mDiffuseBRDF{ std::make_shared<Lambertian>() },
    mAmbientBRDF{ std::make_shared<Lambertian>() }
{}

Matte::Matte(float kd, float ka, Colour color) : Matte{}
{
    setDiffuseReflection(kd);
    setAmbientReflection(ka);
    setDiffuseColour(color);
}

void Matte::setDiffuseReflection(float k)
{
    mDiffuseBRDF->setDiffuseReflection(k);
}

void Matte::setAmbientReflection(float k)
{
    mAmbientBRDF->setDiffuseReflection(k);
}

void Matte::setDiffuseColour(Colour colour)
{
    mDiffuseBRDF->setDiffuseColour(colour);
    mAmbientBRDF->setDiffuseColour(colour);
}

Colour Matte::shade(ShadeRec& sr)
{
    using atlas::math::Ray;
    using atlas::math::Vector;

    Vector wo = -sr.ray.o;
    Colour L = mAmbientBRDF->rho(sr, wo) * sr.world->ambient->L(sr);
    size_t numLights = sr.world->lights.size();

    for (size_t i{ 0 }; i < numLights; ++i)
    {
        Vector wi = sr.world->lights[i]->getDirection(sr);
        float nDotWi = glm::dot(sr.normal, wi);

        if (nDotWi > 0.0f)
        {
            Ray<atlas::math::Vector> shadowRay{};
            shadowRay.o = sr.hitPoint;
            shadowRay.d = glm::normalize(wi);
            bool inShadow = sr.world->lights[i]->inShadow(shadowRay, sr);
            if (!inShadow) {
                L += mDiffuseBRDF->fn(sr, wo, wi) * sr.world->lights[i]->L(sr) *
                    nDotWi;
            }
        }
    }

    return L;
}

Colour Matte::path_shade(ShadeRec& sr, std::shared_ptr<World> world) {
    using atlas::math::Vector;

    Vector wi;
    Vector wo = -sr.ray.d;
    float pdf;
    //Colour f = mDiffuseBRDF->fn
    Colour f = mDiffuseBRDF->sample_f(sr, wo, wi, pdf, world);
    float nDotWi = glm::dot(sr.normal, wi);

    Ray reflect(sr.hitPoint, wi);
    return ((f + world->tracer->trace_ray(reflect, sr.depth + 1, world)) * nDotWi / pdf);
}

// ***** Directional function members *****
Directional::Directional() : Light{}
{}

Directional::Directional(atlas::math::Vector const& d) : Light{}
{
    setDirection(d);
}

void Directional::setDirection(atlas::math::Vector const& d)
{
    mDirection = glm::normalize(d);
}

bool Directional::inShadow([[maybe_unused]] const atlas::math::Ray<atlas::math::Vector> ray, [[maybe_unused]] const ShadeRec& sr) const
{
    using atlas::math::Ray;

    size_t numObjects = sr.world->scene.size();
    float temp = -1;
    for (size_t i{ 0 }; i < numObjects; ++i)
    {
        if (sr.world->scene[i]->shadowHit(ray, temp) && temp > 0) {
            return true;
        }
    }

    return false;
}


atlas::math::Vector Directional::getDirection([[maybe_unused]] ShadeRec& sr)
{
    return mDirection;
}


// ***** Point function members *****
PointL::PointL() : Light{}
{}

PointL::PointL(atlas::math::Vector const& p) : Light{}
{
    setPoint(p);
}

atlas::math::Vector PointL::getDirection([[maybe_unused]] ShadeRec& sr)
{
    return glm::normalize(point - sr.hitPoint);
}

void PointL::setPoint(atlas::math::Point mPoint) {
    point = mPoint;
}

bool PointL::inShadow(const atlas::math::Ray<atlas::math::Vector> ray, const ShadeRec& sr) const {
    using atlas::math::Ray;
    using atlas::math::Point;

    float t{ 10000 };
    size_t num_objects = sr.world->scene.size();
    double d = pow(
        pow((point[0] - ray.o[0]), 2) +
        pow((point[1] - ray.o[1]), 2) +
        pow((point[2] - ray.o[2]), 2), 0.5);

    for (size_t i{ 0 }; i < num_objects; i++) {
        if (sr.world->scene[i]->shadowHit(ray, t) && t < d) {
            return true;
        }
    }
    return false;
}

// ***** Ambient function members *****
Ambient::Ambient() : Light{}
{}

atlas::math::Vector Ambient::getDirection([[maybe_unused]] ShadeRec& sr)
{
    return atlas::math::Vector{ 0.0f };
}

bool Ambient::inShadow([[maybe_unused]] const atlas::math::Ray<atlas::math::Vector> ray, [[maybe_unused]] const ShadeRec& sr) const {
    return false;
}



// ******* Driver Code *******

int main()
{
    using atlas::math::Point;
    using atlas::math::Ray;
    using atlas::math::Vector;

    std::shared_ptr<World> world{ std::make_shared<World>() };

    world->width = 1000;
    world->height = 1000;
    world->background = { 0, 0, 0 };
    world->sampler = std::make_shared<Jittered>(4, 83);

    /*world->scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ -110, 32, -400 }, 100.0f));
    world->scene[0]->setMaterial(
        std::make_shared<Phong>(0.50f, 0.05f, Colour{ 1, 0, 0 }, 15.0f));
    world->scene[0]->setColour({ 1, 1, 1 });*/

    world->scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ 128, 32, -600 }, 100.0f));
    world->scene[0]->setMaterial(
        std::make_shared<Reflective>(0.50f, 0.05f, Colour{ 1, 0, 0 }, 25.0f));
    world->scene[0]->setColour({ 1, 0, 0 });

    world->scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ -20, 32, -400 }, 100.0f));
    world->scene[1]->setMaterial(
        std::make_shared<Transparent>(0.50f, 0.05f, Colour{ 0.01, 0.01, 0.01 }, 2000.0f, Colour{ 1, 1, 1 }, 1.01f));
    world->scene[1]->setColour({ 1, 1, 1 });

    /*world->scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ -128, 32, -700 }, 64.0f));
    world->scene[2]->setMaterial(
        std::make_shared<Matte>(0.60f, 0.2f, Colour{ 0, 1, 0 }));
    world->scene[2]->setColour({ 0, 1, 0 });*/

    atlas::math::Point p1 = { -50, 100, -300 };
    atlas::math::Point p2 = { 100, 150, -150 };
    atlas::math::Point p3 = { 50, 100, -300 };

    world->scene.push_back(
        std::make_shared<Triangle>(p1,
            p2, p3)
    );
    world->scene[2]->setMaterial(
        std::make_shared<Phong>(0.50f, 0.05f, Colour{ 0, 1, 0 }, 35.0f));
    world->scene[2]->setColour({ 0, 1, 1 });

    world->scene.push_back(
        std::make_shared<Plane>(atlas::math::Point{ 0, 150, 0 }, atlas::math::Vector{ 0, 1, 0 }));
    world->scene[3]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 1, 1, 1 }));
    world->scene[3]->setColour({ 0, 0, 1 });


    world->ambient = std::make_shared<Ambient>();
    world->ambient->setColour({ 1, 1, 1 });
    world->ambient->scaleRadiance(1.0f);


    world->lights.push_back(
        std::make_shared<PointL>(PointL({ 150, -150, -450 }))
    );
    world->lights[0]->setColour({ 1,1,1 });
    world->lights[0]->scaleRadiance(6.0f);


    world->lights.push_back(
        std::make_shared<Directional>(Directional{ { -1, 0, 10 } }));
    world->lights[1]->setColour({ 1, 1, 1 });
    world->lights[1]->scaleRadiance(5.0f);

    Point samplePoint{}, pixelPoint{};
    Ray<atlas::math::Vector> ray{ {0, 0, 0}, {0, 0, -1} };

    Pinhole camera{};

    world->tracer = std::make_shared<WhittedTracer>();
    world->maxDepth = 5;

    camera.setEye({ 0.0f, 0.0f, 20.0f });

    camera.computeUVW();

    camera.renderScene(world);

    saveToFile("render.bmp", world->width, world->height, world->image);

    return 0;
}

Colour outOfGamut(Colour pixelAverage) {
    Colour temp = pixelAverage;
    float max = pixelAverage.r;
    if (pixelAverage.g > max) {
        max = pixelAverage.g;
    }
    if (pixelAverage.b > max) {
        max = pixelAverage.b;
    }
    if (max > 1) {
        return { pixelAverage.r / max, pixelAverage.g / max, pixelAverage.b / max };
    }
    else {
        return pixelAverage;
    }
}

void saveToFile(std::string const& filename,
    std::size_t width,
    std::size_t height,
    std::vector<Colour> const& image)
{
    std::vector<unsigned char> data(image.size() * 3);

    for (std::size_t i{ 0 }, k{ 0 }; i < image.size(); ++i, k += 3)
    {
        Colour pixel = image[i];
        data[k + 0] = static_cast<unsigned char>(pixel.r * 255);
        data[k + 1] = static_cast<unsigned char>(pixel.g * 255);
        data[k + 2] = static_cast<unsigned char>(pixel.b * 255);
    }

    stbi_write_bmp(filename.c_str(),
        static_cast<int>(width),
        static_cast<int>(height),
        3,
        data.data());
}